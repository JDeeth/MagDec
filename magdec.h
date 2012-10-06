#ifndef MAGDEC_H
#define MAGDEC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

int my_isnan(double d) {
  return (d != d);              /* IEEE: only NaN is not equal to itself */
}

#define NaN log(-1.0)
const double FT2KM = 1.0/0.0003048;
const double PI = 3.141592654;
const double RAD2DEG = (180.0/PI);

const double EXT_COEFF1 = 0, EXT_COEFF2 = 0, EXT_COEFF3 = 0;

const int MAXDEG = 13;
const int MAXCOEFF = (MAXDEG*(MAXDEG+2)+1); //index starts with 1!, (from old Fortran?)

#ifndef SEEK_SET //for file loading purposes
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

const int IEXT = 0;
const int FALSE = 0;
const int TRUE = 1;
const int RECL = 81; //expected size of line in model (including trailing \0 evidently)

//#define MAXINBUFF RECL+14
const int MAXINBUFF = RECL + 14;
/** Max size of in buffer **/

//#define MAXREAD MAXINBUFF-2
const int MAXREAD = MAXINBUFF - 2;
/** Max to read 2 less than total size (just to be safe) **/

//#define MAXMOD 30
const int MAXMOD = 30;
/** Max number of models in a file **/

//#define PATH MAXREAD
const int PATH = MAXREAD;
/** Max path and filename length **/

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
*                           Subroutine getshc                              *
*                                                                          *
****************************************************************************
*                                                                          *
*     Reads spherical harmonic coefficients from the specified             *
*     model into an array.                                                 *
*                                                                          *
*     Input:                                                               *
*           stream     - Logical unit number                               *
*           iflag      - Flag for SV equal to ) or not equal to 0          *
*                        for designated read statements                    *
*           strec      - Starting record number to read from model         *
*           nmax_of_gh - Maximum degree and order of model                 *
*                                                                          *
*     Output:                                                              *
*           gh1 or 2   - Schmidt quasi-normal internal spherical           *
*                        harmonic coefficients                             *
*                                                                          *
*     Return values:                                                       *
*           -2: first two columns of input.cof not in correct pattern      *
*           0: success, or no input stream                                 *
*                                                                          *
*     FORTRAN                                                              *
*           Bill Flanagan                                                  *
*           NOAA CORPS, DESDIS, NGDC, 325 Broadway, Boulder CO.  80301     *
*                                                                          *
*     C                                                                    *
*           C. H. Shaffer                                                  *
*           Lockheed Missiles and Space Company, Sunnyvale CA              *
*           August 15, 1988                                                *
*                                                                          *
***************************************************************************/

int getshc(FILE *stream, char file[PATH], int iflag, long strec, int nmax_of_gh, double* gh) {
  char  inbuff[MAXINBUFF];
  int ii,m,n,mm,nn;
  double g,hh;

  stream = fopen(file, "rt");
  if (stream == NULL) {
    printf("\nError on opening file %s", file);
  }
  else {
    ii = 0;
    fseek(stream,strec,SEEK_SET); //find beginning of model in input.cof

    for ( nn = 1; nn <= nmax_of_gh; ++nn) {
      for (mm = 0; mm <= nn; ++mm) {
        if (iflag == 1) { //read first pair of g,h columns
          fgets(inbuff, MAXREAD, stream);
          sscanf(inbuff, "%d %d %lg %lg %*lg %*lg %*s %*d",
                 &n, &m, &g, &hh); //%*x is identified and ignored
        }
        else { //read second pair of g,h columns
          fgets(inbuff, MAXREAD, stream);
          sscanf(inbuff, "%d %d %*lg %*lg %lg %lg %*s %*d",
                 &n, &m, &g, &hh);
        }

        if ((nn != n) || (mm != m)) { //this serves to check n,m are increasing in correct
          fclose(stream);             //sequence in input.cof
          return(-2);
        }

        ++ii;
        gh[ii] = g;

        if (m != 0) {
          ++ii;
          gh[ii] = hh;
        }
      }
    }
  }
  fclose(stream);
  return(0);
}

/***************************************************************************
*                                                                          *
*                           Subroutine extrapsh                            *
*                                                                          *
****************************************************************************
*                                                                          *
*     Extrapolates linearly a spherical harmonic model with a              *
*     rate-of-change model.                                                *
*                                                                          *
*     Input:                                                               *
*           date     - date of resulting model (in decimal year)           *
*           dte1     - date of base model                                  *
*           nmax1    - maximum degree and order of base model              *
*           gh1      - Schmidt quasi-normal internal spherical             *
*                      harmonic coefficients of base model                 *
*           nmax2    - maximum degree and order of rate-of-change model    *
*           gh2      - Schmidt quasi-normal internal spherical             *
*                      harmonic coefficients of rate-of-change model       *
*                                                                          *
*     Output:                                                              *
*           gha or b - Schmidt quasi-normal internal spherical             *
*                    harmonic coefficients                                 *
*           nmax   - maximum degree and order of resulting model           *
*                                                                          *
*     FORTRAN                                                              *
*           A. Zunde                                                       *
*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     *
*                                                                          *
*     C                                                                    *
*           C. H. Shaffer                                                  *
*           Lockheed Missiles and Space Company, Sunnyvale CA              *
*           August 16, 1988                                                *
*                                                                          *
***************************************************************************/

int extrapsh(const double &date, const double &dte1, const int &nmax1, const int &nmax2,
             double* gh, const double*const gh1, const double*const gh2)
{
  int   nmax;
  int   k, l;
  int   ii;
  double factor; // decimal years since epoch of base model

  factor = date - dte1;

  if (nmax1 == nmax2) {
    k =  nmax1 * (nmax1 + 2);
    nmax = nmax1;
  }

  else {
    if (nmax1 > nmax2) {
      k = nmax2 * (nmax2 + 2);
      l = nmax1 * (nmax1 + 2);
      for ( ii = k + 1; ii <= l; ++ii) {
        gh[ii] = gh1[ii];
      }
      nmax = nmax1;
    }
    else {
      k = nmax1 * (nmax1 + 2);
      l = nmax2 * (nmax2 + 2);
      for ( ii = k + 1; ii <= l; ++ii) {
        gh[ii] = factor * gh2[ii];
      }
      nmax = nmax2;
    }
  }

  for ( ii = 1; ii <= k; ++ii) {
    gh[ii] = gh1[ii] + factor * gh2[ii];
  }

  return(nmax);
}

/***************************************************************************
*                                                                          *
*                           Subroutine interpsh                            *
*                                                                          *
****************************************************************************
*                                                                          *
*     Interpolates linearly, in time, between two spherical harmonic       *
*     models.                                                              *
*                                                                          *
*     Input:                                                               *
*           date     - date of resulting model (in decimal year)           *
*           dte1     - date of earlier model                               *
*           nmax1    - maximum degree and order of earlier model           *
*           gh1      - Schmidt quasi-normal internal spherical             *
*                      harmonic coefficients of earlier model              *
*           dte2     - date of later model                                 *
*           nmax2    - maximum degree and order of later model             *
*           gh2      - Schmidt quasi-normal internal spherical             *
*                      harmonic coefficients of internal model             *
*                                                                          *
*     Output:                                                              *
*           gha or b - coefficients of resulting model                     *
*           nmax     - maximum degree and order of resulting model         *
*                                                                          *
*     FORTRAN                                                              *
*           A. Zunde                                                       *
*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     *
*                                                                          *
*     C                                                                    *
*           C. H. Shaffer                                                  *
*           Lockheed Missiles and Space Company, Sunnyvale CA              *
*           August 17, 1988                                                *
*                                                                          *
***************************************************************************/


int interpsh(const double &date, const double &dte1, const int &nmax1,
             const double &dte2, const int &nmax2,
             double* gh, const double*const gh1, const double*const gh2)
{
  int   nmax;
  int   k, l;
  int   ii;
  double factor;

  factor = (date - dte1) / (dte2 - dte1);

  if (nmax1 == nmax2) {
    k =  nmax1 * (nmax1 + 2);
    nmax = nmax1;
  }

  else {
    if (nmax1 > nmax2) {
      k = nmax2 * (nmax2 + 2);
      l = nmax1 * (nmax1 + 2);

      for ( ii = k + 1; ii <= l; ++ii) {
        gh[ii] = gh1[ii] + factor * (-gh1[ii]);
      }

      nmax = nmax1;
    }
    else {
      k = nmax1 * (nmax1 + 2);
      l = nmax2 * (nmax2 + 2);

      for ( ii = k + 1; ii <= l; ++ii) {
        gh[ii] = factor * gh2[ii];
      }

      nmax = nmax2;
    }
  }

  for ( ii = 1; ii <= k; ++ii) {
    gh[ii] = gh1[ii] + factor * (gh2[ii] - gh1[ii]);
  }

  return(nmax);
}

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
  int npq;          //number of paired qoefficients? (gh pairs)

  pointComponents xyz; //for output
  xyz.x = 0;
  xyz.y = 0;
  xyz.z = 0;

  a2 = 40680631.59;            /* WGS84 */
  b2 = 40408299.98;            /* WGS84 */
  r = point.elev;
  slat = sin( point.lat * dtr );

  if ((90.0 - point.lat) < 0.001) {
    aa = 89.999;            /*  300 ft. from North pole  */
  }
  else {
    if ((90.0 + point.lat) < 0.001) {
      aa = -89.999;        /*  300 ft. from South pole  */
    }
    else {
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

  for ( k = 1; k <= npq; ++k) { //k is iterating through g,h pairs

    if (n < m) {
      m = 0;
      n++;
      rr = pow(ratio, n + 2);
      fn = n;
    }

    fm = m;

    if (k >= 5) {
      if (m == n) {
        aa = sqrt( 1.0 - 0.5/fm );
        j = k - n - 1;
        p[k] = (1.0 + 1.0/fm) * aa * clat * p[j];
        q[k] = aa * (clat * q[j] + slat/fm * p[j]);
        sl[m] = sl[m-1] * cl[1] + cl[m-1] * sl[1];
        cl[m] = cl[m-1] * cl[1] - sl[m-1] * sl[1];
      }
      else {
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

    if (m == 0) {
      xyz.x = xyz.x + aa * q[k];
      xyz.z = xyz.z - aa * p[k];
      l++;
    }
    else {
      bb = rr * ghArray[l+1];
      cc = aa * cl[m] + bb * sl[m];
      xyz.x = xyz.x + cc * q[k];
      xyz.z = xyz.z - cc * p[k];

      if (clat > 0) {
        xyz.y = xyz.y + (aa * sl[m] - bb * cl[m]) *
                fm * p[k]/((fn + 1.0) * clat);
      }
      else {
        xyz.y = xyz.y + (aa * sl[m] - bb * cl[m]) * q[k] * slat;
      }

      l += 2;
    }

    m++;
  }

  if (iext != 0) {
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

/***************************************************************************
*                                                                          *
*                           Subroutine dihf                                *
*                                                                          *
****************************************************************************
*                                                                          *
*     Computes the geomagnetic d, i, h, and f from x, y, and z.            *
*                                                                          *
*     Input:                                                               *
*           x  - northward component                                       *
*           y  - eastward component                                        *
*           z  - vertically-downward component                             *
*                                                                          *
*     Output:                                                              *
*           d  - declination                                               *
*           i  - inclination                                               *
*           h  - horizontal intensity                                      *
*           f  - total intensity                                           *
*                                                                          *
*     FORTRAN                                                              *
*           A. Zunde                                                       *
*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     *
*                                                                          *
*     C                                                                    *
*           C. H. Shaffer                                                  *
*           Lockheed Missiles and Space Company, Sunnyvale CA              *
*           August 22, 1988                                                *
*                                                                          *
***************************************************************************/

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
