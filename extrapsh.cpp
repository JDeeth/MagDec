#include "magdec.h"

int main () {
  double gha[MAXCOEFF];
  for(int i = 1; i<=MAXCOEFF; ++i) {
      gha[i]=i;
    }
  extrapsh(0,0,0,0,&gha[0],&gha[0],&gha[0]);
  return 0;
}
