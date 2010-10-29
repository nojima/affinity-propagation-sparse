// An Implementation of Affinity Propergation
// See: Clustering by Passing Messages Between Data Points

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>
#include "ap.h"
using namespace std;

namespace {
  struct Edge {
    int src;      // index of source
    int dst;      // index of destination
    double s;     // similarity s(src, dst)
    double r;     // responsibility r(src, dst)
    double a;     // availability a(src, dst)

    Edge(int src, int dst, double s, double r, double a)
      : src(src), dst(dst), s(s), r(r), a(a) {}
  };

  typedef vector<Edge*> Edges;

  struct Graph {
    int n;                // the number of vertices
    Edges* outEdges;      // array of out edges of corresponding vertices
    Edges* inEdges;       // array of in edges of corresponding vertices
  };

  // Build graph from sparse similarity matrix stored in COO format.
  // Input specification is following:
  // First line contains an integer standing for the size of the matrix.
  // Following lines each contain two integers and a real number standing for
  // the row index i, the column index j and the similarity s(i,j) respectively.
  // Input ends with an end-of-file.
  // Note that this function does not check any errors in the given input.
  // Parameter:
  //   input: Input file handle.
  //   prefType:
  //     1: use median of similarities as preference
  //     2: use minimum of similarities as preference
  //     3: use min - (max - min) of similarities as preference
  Graph* buildGraph(FILE* input, int prefType)
  {
    Graph* graph = new Graph;
    fscanf(input, "%d", &graph->n);
    graph->outEdges = new Edges[graph->n];
    graph->inEdges = new Edges[graph->n];

    // read similarity matrix
    int i, j;
    double s;
    vector<double> ss;
    while (fscanf(input, "%d%d%lf", &i, &j, &s) != EOF) {
      if (i == j) { continue; }
      // add small noise to avoid degeneracies
      s += (1e-16 * s + 1e-300) * (rand() / (RAND_MAX + 1.0));
      // make new edge
      Edge* e = new Edge(i, j, s, 0, 0);
      graph->outEdges[i].push_back(e);
      graph->inEdges[j].push_back(e);
      // store similarities for preference calculation
      ss.push_back(s);
    }

    // calculate preferences
    double pref;
    if (prefType == 1) {
      sort(ss.begin(), ss.end());
      int m = ss.size();
      pref = (m % 2) ? ss[m/2] : (ss[m/2 - 1] + ss[m/2]) / 2.0;
    } else if (prefType == 2) {
      pref = *min_element(ss.begin(), ss.end());
    } else if (prefType == 3) {
      double minValue = *min_element(ss.begin(), ss.end());
      double maxValue = *max_element(ss.begin(), ss.end());
      pref = 2*minValue - maxValue;
    } else {
      assert(false);      // invalid prefType
    }
    for (int i = 0; i < graph->n; ++i) {
      Edge* e = new Edge(i, i, pref, 0, 0);
      graph->outEdges[i].push_back(e);
      graph->inEdges[i].push_back(e);
    }

    return graph;
  }

  void destroyGraph(Graph* g)
  {
    for (int i = 0; i < g->n; ++i) {
      for (size_t j = 0; j < g->outEdges[i].size(); ++j) {
        delete g->outEdges[i][j];
      }
    }
    delete [] g->outEdges;
    delete [] g->inEdges;
    delete g;
  }

  inline void update(double& variable, double newValue, double damping)
  {
    variable = damping * variable + (1.0 - damping) * newValue;
  }

  void updateResponsibilities(Graph* graph, double damping)
  {
    for (int i = 0; i < graph->n; ++i) {
      Edges& edges = graph->outEdges[i];
      int m = edges.size();
      double max1 = -HUGE_VAL, max2 = -HUGE_VAL;
      double argmax1 = -1;
      for (int k = 0; k < m; ++k) {
        double v = edges[k]->s + edges[k]->a;
        if (v > max1) { swap(max1, v); argmax1 = k; }
        if (v > max2) { max2 = v; }
      }
      // update responsibilities
      for (int k = 0; k < m; ++k) {
        if (k != argmax1) {
          update(edges[k]->r, edges[k]->s - max1, damping);
        } else {
          update(edges[k]->r, edges[k]->s - max2, damping);
        }
      }
    }
  }

  void updateAvailabilities(Graph* graph, double damping)
  {
    for (int k = 0; k < graph->n; ++k) {
      Edges& edges = graph->inEdges[k];
      int m = edges.size();
      // calculate sum of positive responsibilities
      double sum = 0.0;
      for (int i = 0; i < m-1; ++i) {
        sum += max(0.0, edges[i]->r);
      }
      // calculate availabilities
      double rkk = edges[m-1]->r;
      for (int i = 0; i < m-1; ++i) {
        update(edges[i]->a, min(0.0, rkk + sum - max(0.0, edges[i]->r)), damping);
      }
      // calculate self-availability
      update(edges[m-1]->a, sum, damping);
    }
  }

  bool updateExamplars(Graph* graph, vector<int>& examplar)
  {
    bool changed = false;
    for (int i = 0; i < graph->n; ++i) {
      Edges& edges = graph->outEdges[i];
      int m = edges.size();
      double maxValue = -HUGE_VAL;
      int argmax = i;
      for (int k = 0; k < m; ++k) {
        double value = edges[k]->a + edges[k]->r;
        if (value > maxValue) {
          maxValue = value;
          argmax = edges[k]->dst;
        }
      }
      if (examplar[i] != argmax) {
        examplar[i] = argmax;
        changed = true;
      }
    }
    return changed;
  }
}

// Cluster data points with Affinity Propagation.
// Parameters:
//   input: Input file which contains sparse similarity matrix. see buildGraph().
//   prefType: Specify what kind of preference we use. see buildGraph().
//   damping: The damping factor. (0.5 <= damping < 1.0)
//   maxit: The maximum number of iterations.
//   convit: Specify how many iterations this algorithm stops when examplars
//           did not change for.
// Returns:
//   Array of examplars of corresponding data points.
vector<int> affinityPropagation(FILE* input, int prefType, double damping, int maxit, int convit)
{
  assert(0.499 < damping && damping < 1.0);

  Graph* graph = buildGraph(input, prefType);
  vector<int> examplar(graph->n, -1);

  for (int i = 0, nochange = 0; i < maxit && nochange < convit; ++i, ++nochange) {
    updateResponsibilities(graph, damping);
    updateAvailabilities(graph, damping);
    if (updateExamplars(graph, examplar)) { nochange = 0; }
  }
  
  destroyGraph(graph);
  return examplar;
}
