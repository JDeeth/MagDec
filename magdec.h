#ifndef MAGDEC_H
#define MAGDEC_H

#include <math.h> 

int my_isnan(double d)
{
  return (d != d);              /* IEEE: only NaN is not equal to itself */
}

#define NaN log(-1.0)
const double FT2KM = 1.0/0.0003048;
const double PI = 3.141592654;
const double RAD2DEG = (180.0/PI);

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


#endif
