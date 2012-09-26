#include "magdec.h"
#include <iostream>

/****************************************************************************/
/*                                                                          */
/*                           Subroutine dihf                                */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Computes the geomagnetic d, i, h, and f from x, y, and z.            */
/*                                                                          */
/*     Input:                                                               */
/*           x  - northward component                                       */
/*           y  - eastward component                                        */
/*           z  - vertically-downward component                             */
/*                                                                          */
/*     Output:                                                              */
/*           d  - declination                                               */
/*           i  - inclination                                               */
/*           h  - horizontal intensity                                      */
/*           f  - total intensity                                           */
/*                                                                          */
/*     FORTRAN                                                              */
/*           A. Zunde                                                       */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 22, 1988                                                */
/*                                                                          */
/****************************************************************************/

pointField dihf (const pointComponents &xyz)
{
  double sn = 0.0001; //small number
  double h2;
  double hpx;

  double x = xyz.x;
  double y = xyz.y;
  double z = xyz.z;

  double d, i, h, f;
  
  h2 = x*x + y*y;
  h = sqrt(h2);       /* calculate horizontal intensity */
  f = sqrt(h2 + z*z);      /* calculate total intensity */
  if (f < sn)
    {
      d = NaN;        /* If d and i cannot be determined, */
      i = NaN;        /*       set equal to NaN         */
    }
  else
    {
      i = atan2(z, h);
      if (h < sn)
        {
          d = NaN;
        }
      else
        {
          hpx = h + x;
          if (hpx < sn) //
            {
              d = PI;
            }
          else
            {
              d = 2.0 * atan2(y, hpx);
            }
        }
    }

  pointField dihf;
  dihf.d = d * RAD2DEG;
  dihf.i = i * RAD2DEG;
  dihf.h = h;
  dihf.f = f;
  return(dihf);
}

//dummy main
using namespace std;
int main () {
  pointComponents input;
  input.x = 27547.2; //field components for 0N, 0E, 0ft, 2012.0
  input.y = -2817.5;
  input.z = -15628.3;
  cout <<   "Input x, y, z:" << endl <<
            input.x << " " <<
            input.y << " " <<
            input.z << endl;

  pointField output = dihf(input);

  cout << "Output d, i, h, f:" << endl <<
          output.d << " " << output.i<< " " <<
          output.h << " " << output.f << endl;

  return 0;
}
