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

const double EXT_COEFF1 = 0, EXT_COEFF2 = 0, EXT_COEFF3 = 0;

const int MAXDEG = 13;
const int MAXCOEFF = (MAXDEG*(MAXDEG+2)+1); //index starts with 1!, (from old Fortran?)

struct pointCoords {
  double lat, lon, elev; //latitude, longitude (degrees),
                         //elevation from WGS84 (km),
                         //in geodetic coordinates
};

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

/***************************************************************************
*                                                                          *
*                           Subroutine shval3                              *
*                                                                          *
****************************************************************************
*                                                                          *
*     Calculates field components from spherical harmonic (sh)             *
*     models.                                                              *
*                                                                          *
*     Input:                                                               *
*           igdgc     - indicates coordinate system used; set equal        *
*                       to 1 if geodetic, 2 if geocentric                  *
*           latitude  - north latitude, in degrees                         *
*           longitude - east longitude, in degrees                         *
*           elev      - WGS84 altitude above ellipsoid (igdgc=1), or       *
*                       radial distance from earth's center (igdgc=2)      *
*           a2,b2     - squares of semi-major and semi-minor axes of       *
*                       the reference spheroid used for transforming       *
*                       between geodetic and geocentric coordinates        *
*                       or components                                      *
*           nmax      - maximum degree and order of coefficients           *
*           iext      - external coefficients flag (=0 if none)            *
*           ext1,2,3  - the three 1st-degree external coefficients         *
*                       (not used if iext = 0)                             *
*                                                                          *
*     Output:                                                              *
*           x         - northward component                                *
*           y         - eastward component                                 *
*           z         - vertically-downward component                      *
*                                                                          *
*     based on subroutine 'igrf' by D. R. Barraclough and S. R. C. Malin,  *
*     report no. 71/1, institute of geological sciences, U.K.              *
*                                                                          *
*     FORTRAN                                                              *
*           Norman W. Peddie                                               *
*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     *
*                                                                          *
*     C                                                                    *
*           C. H. Shaffer                                                  *
*           Lockheed Missiles and Space Company, Sunnyvale CA              *
*           August 17, 1988                                                *
*                                                                          *
***************************************************************************/

pointComponents shval3(
    pointCoords point, //geodetic coordinates
    int nmax, //"Maximum degree and order of coefficients
    const double* ghArray,
    int iext, double ext1, double ext2, double ext3)
{
  double earths_radius = 6371.2; //in km
  double dtr = 0.01745329;
  double slat;
  double clat;
  double ratio;
  double aa, bb, cc, dd;
  double sd;
  double cd;
  double r;
  double a2;
  double b2;
  double rr;
  double fm,fn;
  double sl[14];
  double cl[14];
  double p[119];
  double q[119];
  int ii,j,k,l,m,n;
  int npq;

  pointComponents xyz; //for output
  xyz.x = 0;
  xyz.y = 0;
  xyz.z = 0;

  a2 = 40680631.59;            /* WGS84 */
  b2 = 40408299.98;            /* WGS84 */
  r = point.elev;
  slat = sin( point.lat * dtr );
  if ((90.0 - point.lat) < 0.001)
    {
      aa = 89.999;            /*  300 ft. from North pole  */
    }
  else
    {
      if ((90.0 + point.lat) < 0.001)
        {
          aa = -89.999;        /*  300 ft. from South pole  */
        }
      else
        {
          aa = point.lat;
        }
    }

  clat = cos( aa * dtr );
  sl[1] = sin( point.lon * dtr );
  cl[1] = cos( point.lon * dtr );

  sd = 0.0;
  cd = 1.0;
  l = 1;
  n = 0;
  m = 1;
  npq = (nmax * (nmax + 3)) / 2;

  aa = a2 * clat * clat;
  bb = b2 * slat * slat;
  cc = aa + bb;
  dd = sqrt( cc );
  r = sqrt( point.elev * (point.elev + 2.0 * dd) + (a2 * aa + b2 * bb) / cc );
  cd = (point.elev + dd) / r;
  sd = (a2 - b2) / dd * slat * clat / r;
  aa = slat;
  slat = slat * cd - clat * sd;
  clat = clat * cd + aa * sd;

  ratio = earths_radius / r;
  aa = sqrt( 3.0 );
  p[1] = 2.0 * slat;
  p[2] = 2.0 * clat;
  p[3] = 4.5 * slat * slat - 1.5;
  p[4] = 3.0 * aa * clat * slat;
  q[1] = -clat;
  q[2] = slat;
  q[3] = -3.0 * clat * slat;
  q[4] = aa * (slat * slat - clat * clat);
  for ( k = 1; k <= npq; ++k)
    {
      if (n < m)
        {
          m = 0;
          n++;
          rr = pow(ratio, n + 2);
          fn = n;
        }
      fm = m;
      if (k >= 5)
        {
          if (m == n)
            {
              aa = sqrt( 1.0 - 0.5/fm );
              j = k - n - 1;
              p[k] = (1.0 + 1.0/fm) * aa * clat * p[j];
              q[k] = aa * (clat * q[j] + slat/fm * p[j]);
              sl[m] = sl[m-1] * cl[1] + cl[m-1] * sl[1];
              cl[m] = cl[m-1] * cl[1] - sl[m-1] * sl[1];
            }
          else
            {
              aa = sqrt( fn*fn - fm*fm);
              bb = sqrt( ((fn - 1.0)*(fn-1.0)) - (fm * fm) )/aa;
              cc = (2.0 * fn - 1.0)/aa;
              ii = k - n;
              j = k - 2 * n + 1;
              p[k] = (fn + 1.0) * (cc * slat/fn * p[ii] - bb/(fn - 1.0) * p[j]);
              q[k] = cc * (slat * q[ii] - clat/fn * p[ii]) - bb * q[j];
            }
        }
      aa = rr * ghArray[l];
      if (m == 0)
        {
          xyz.x = xyz.x + aa * q[k];
          xyz.z = xyz.z - aa * p[k];
          l++;
        }
      else
        {
          bb = rr * ghArray[l+1];
          cc = aa * cl[m] + bb * sl[m];
          xyz.x = xyz.x + cc * q[k];
          xyz.z = xyz.z - cc * p[k];
          if (clat > 0)
            {
              xyz.y = xyz.y + (aa * sl[m] - bb * cl[m]) *
                  fm * p[k]/((fn + 1.0) * clat);
            }
          else
            {
              xyz.y = xyz.y + (aa * sl[m] - bb * cl[m]) * q[k] * slat;
            }
          l = l + 2;
        }
      m++;
    }
  if (iext != 0)
    {
      aa = ext2 * cl[1] + ext3 * sl[1];
      xyz.x = xyz.x - ext1 * clat + aa * slat;
      xyz.y = xyz.y + ext2 * sl[1] - ext3 * cl[1];
      xyz.z = xyz.z + ext1 * slat + aa * clat;
    }

  aa = xyz.x;
  xyz.x = xyz.x * cd + xyz.z * sd;
  xyz.z = xyz.z * cd - aa * sd;

  return(xyz);
}

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
