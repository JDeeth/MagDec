#include "magdec.h"
#include <iostream>

using namespace std;

int main() {
  FILE *stream;
  char file[PATH] = "igrf11.cof\0";
  int iflag = 5;
  long strec = 2;
  int nmax_of_gh = 2;
  double gh[255];
  getshc(stream, &file[0], iflag, strec, nmax_of_gh, &gh[0]);
  for (int i=0; i<20; i++) {
      cout << endl << gh[i];
    }
  return 0;
}
