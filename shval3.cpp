#include "magdec.h"

int main () {
  double inArray[MAXCOEFF];
  for (int i=1; i<=MAXCOEFF; i++) {
      inArray[i] = i;
    }
  shval3(
      0,0,0,
      10,
      &inArray[0],
      3, 10.0, 10.0, 10.0);

  return 0;
}
