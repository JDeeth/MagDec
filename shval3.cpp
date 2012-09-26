#include "magdec.h"

int main () {
  double inArray[MAXCOEFF];
  for (int i=1; i<=MAXCOEFF; i++) {
      inArray[i] = i;
    }
  pointCoords loc;
  loc.lat = 0;
  loc.lon = 0;
  loc.elev = 0;
  shval3(
      loc,
      10,
      &inArray[0],
      3, 10.0, 10.0, 10.0);

  return 0;
}
