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

struct pointComponents {
  double x, y, z; //north, east, and downward component of magnetic field at
  //a [geocentric/geodetic?] point. Units are nT.
};

struct pointField {
  double d; //declination from geographic north (deg)
  double i; //inclination (deg)
  double h; //horizontal field strength
  double f; //total field strength
};

pointField dihf (const pointComponents &xyz)
{
  int j;
  double sn; //small number
  double h2;
  double hpx;

  double x = xyz.x;
  double y = xyz.y;
  double z = xyz.z;

  double d, i, h, f;

  sn = 0.0001;
  
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
  dihf.d = d;
  dihf.i = i;
  dihf.h = h;
  dihf.f = f;
  return(dihf);
}

//dummy main
using namespace std;
int main () {
  pointComponents input;
  input.x = 27547.2;
  input.y = -2817.5;
  input.z = -15628.3;
  cout <<   "Input x, y, z:" << endl <<
            input.x << " " <<
            input.y << " " <<
            input.z << endl;

  pointField output = dihf(input);

  //copypasta
  int   ddeg,ideg;
  double dmin,imin;
  double d, i;
  d = output.d;
  i = output.i;

  /* Change d and i to deg and min */
  ddeg=(int)d;
  dmin=(d-(double)ddeg)*60;
  if (d > 0 && dmin >= 59.5)
    {
      dmin -= 60.0;
      ddeg++;
    }
  if (d < 0 && dmin <= -59.5)
    {
      dmin += 60.0;
      ddeg--;
    }

  if (ddeg!=0) dmin=fabs(dmin);
  ideg=(int)i;
  imin=(i-(double)ideg)*60;
  if (i > 0 && imin >= 59.5)
    {
      imin -= 60.0;
      ideg++;
    }
  if (i < 0 && imin <= -59.5)
    {
      imin += 60.0;
      ideg--;
    }
  if (ideg!=0) imin=fabs(imin);

  //copypasta end

  cout << "Output d, i, h, f:" << endl <<
          ddeg << "\'" << dmin << " " <<
          ideg << "\'" << imin << " " <<
          output.h << " " << output.f << endl;

  return 0;
}
