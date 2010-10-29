#include <cstdio>
#include <vector>
#include "ap.h"
using namespace std;

// test
int main()
{
  vector<int> examplar = affinityPropagation(stdin);
  printf("examplar:");
  for (size_t i = 0; i < examplar.size(); ++i) {
    printf(" %d", examplar[i]);
  }
  puts("");
  return 0;
}
