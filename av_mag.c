/****************************************************************************/
/*                                                                          */
/*       Disclaimer: This C version of the Geomagnetic Field                */
/*       Modeling software is being supplied to aid programmers in          */
/*       integrating spherical harmonic field modeling calculations         */
/*       into their code. It is being distributed unoffically. The          */
/*       National Geophysical Data Center does not support it as            */
/*       a data product or guarantee it's correctness. The limited          */
/*       testing done on this code seems to indicate that the results       */
/*       obtained from this program are compariable to those obtained       */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     In regards to the disclaimer, to the best of my knowlege and quick   */
/*     testing, this program, generates numaric output which is withing .1  */
/*     degress of the current fortran version.  However, it *is* a program  */
/*     and most likely contains bugs.                                       */
/*                                              dio  6-6-96                 */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     This is version 3.01 of the source code but still represents the     */
/*     geomag30 executable.                                                 */
/*                                              dio 9-17-96                 */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Bug fix- the range of dates error.  There was a difference in the    */
/*     output between the last value on a range of dates and the individual */
/*     value for the same date.  Lets make this version 3.02                */
/*                                              ljh 11-20-98                */
/*                                                                          */
/****************************************************************************/
/*      Program was modified so that it can accept year 2000 models         */
/*      This required that the number of blank spaces on  model header      */
/*      records be decreased from four to three                             */
/*                                                                          */
/*      This program calculates the geomagnetic field values from           */
/*      a spherical harmonic model.  Inputs required by the user are:       */
/*      a spherical harmonic model data file, coordinate preference,        */
/*      elevation, date/range-step, latitude, and longitude.                */
/*                                                                          */
/*         Spherical Harmonic                                               */
/*         Model Data File       :  Name of the data file containing the    */
/*                                  spherical harmonic coefficients of      */
/*                                  the chosen model.  The model and path   */
/*                                  must be less than PATH chars.           */
/*                                                                          */
/*         Coordinate Preference :  Geodetic (measured from                 */
/*                                  the surface of the earth),              */
/*                                  or geocentric (measured from the        */
/*                                  center of the earth).                   */
/*                                                                          */
/*         Elevation             :  Elevation above sea level in kilometers.*/
/*                                  if geocentric coordinate preference is  */
/*                                  used then the elevation must be in the  */
/*                                  range of 6370.20 km - 6971.20 km as     */
/*                                  measured from the center of the earth.  */
/*                                  Enter elevation in kilometers in the    */
/*                                  form xxxx.xx.                           */
/*                                                                          */
/*         Date                  :  Date, in decimal years, for which to    */
/*                                  calculate the values of the magnetic    */
/*                                  field.  The date must be within the     */
/*                                  limits of the model chosen.             */
/*                                                                          */
/*         Latitude              :  Entered in decimal degrees in the       */
/*                                  form xxx.xxx.  Positive for northern    */
/*                                  hemisphere, negative for the southern   */
/*                                  hemisphere.                             */
/*                                                                          */
/*         Longitude             :  Entered in decimal degrees in the       */
/*                                  form xxx.xxx.  Positive for eastern     */
/*                                  hemisphere, negative for the western    */
/*                                  hemisphere.                             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*      Subroutines called :  degrees_to_decimal,julday,getshc,interpsh,    */
/*                            extrapsh,shval3,dihf,safegets                 */
/*                                                                          */
/****************************************************************************/

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "av_mag.h"
#include "av_geo.h"

#define MAXDEG (13)
#define MAXCOEFF (MAXDEG*(MAXDEG+2))
struct GeoContext {
  float gh1[MAXCOEFF];
  float gh2[MAXCOEFF];
  float gha[MAXCOEFF];                       /* Geomag global variables */
  float ghb[MAXCOEFF];
  float d, f, h, i;
  float dtemp,ftemp,htemp,itemp;
  float x, y, z;
  float xtemp,ytemp,ztemp;
};

void
geomag(struct GeoContext *pctx, const struct GeoModel *pmodel, float sdate,
       float latitude, float longitude, float elev);


static void
getshc(int iflag, const struct GeoData *pdata, int nmax_of_gh, float *pgh);

static int
extrapsh(struct GeoContext *pctx, float date, float dte1, int nmax1,
	 int nmax2, float *pgh);

static int
interpsh(struct GeoContext *pctx, float date, float dte1, int nmax1,
	 float dte2, int nmax2, float *pgh);

static int
shval3(struct GeoContext *pctx, float flat, float flon, float elev, int nmax, int gh);

static void
dihf3(struct GeoContext *pctx);

static void
dihf4(struct GeoContext *pctx);

float
magnetic_variation(float sdate, float latitude,
		   float longitude, float elev) {
  struct GeoContext ctx;
  unsigned i;

  for (i = 1; i < GEO_MODEL_COUNT; i++)
    if (sdate < GEO_MODELS[i].yrmin)
      break;

  geomag(&ctx, GEO_MODELS + i - 1, sdate, latitude, longitude, elev);

  return -ctx.d;
}

void
geomag(struct GeoContext *pctx, const struct GeoModel *pmodel, float sdate,
       float latitude, float longitude, float elev) {
  int nmax;

  /** This will compute everything needed for 1 point in time. **/

  if(pmodel->max2 == 0) {
    getshc(1, pmodel->data, pmodel->max1, pctx->gh1);
    getshc(1, (pmodel + 1)->data, (pmodel + 1)->max1, pctx->gh2);
    nmax = interpsh(pctx, sdate, pmodel->yrmin, pmodel->max1,
		    (pmodel + 1)->yrmin, (pmodel + 1)->max1, pctx->gha);
    nmax = interpsh(pctx, sdate+1, pmodel->yrmin , pmodel->max1,
		    (pmodel + 1)->yrmin, (pmodel + 1)->max1, pctx->ghb);
  } else {
    getshc(1, pmodel->data, pmodel->max1, pctx->gh1);
    getshc(0, pmodel->data, pmodel->max2, pctx->gh2);
    nmax = extrapsh(pctx, sdate, pmodel->yrmin, pmodel->max1, pmodel->max2, pctx->gha);
    nmax = extrapsh(pctx, sdate+1, pmodel->yrmin, pmodel->max1, pmodel->max2, pctx->ghb);
  }

  /* Do the first calculations */
  shval3(pctx, latitude, longitude, elev / 1000.0f, nmax, 3);

  dihf3(pctx);
  shval3(pctx, latitude, longitude, elev / 1000.0f, nmax, 4);
  dihf4(pctx);

  /** Above will compute everything for 1 point in time.  **/
}

/****************************************************************************/
/*                                                                          */
/*                           Subroutine getshc                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Reads spherical harmonic coefficients from the specified             */
/*     model into an array.                                                 */
/*                                                                          */
/*     Input:                                                               */
/*           iflag      - Flag for SV equal to ) or not equal to 0          */
/*                        for designated read statements                    */
/*           strec      - Starting record number to read from model         */
/*           nmax_of_gh - Maximum degree and order of model                 */
/*                                                                          */
/*     Output:                                                              */
/*           gh1 or 2   - Schmidt quasi-normal internal spherical           */
/*                        harmonic coefficients                             */
/*                                                                          */
/*     FORTRAN                                                              */
/*           Bill Flanagan                                                  */
/*           NOAA CORPS, DESDIS, NGDC, 325 Broadway, Boulder CO.  80301     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 15, 1988                                                */
/*                                                                          */
/****************************************************************************/


static void
getshc(int iflag, const struct GeoData *pdata, int nmax_of_gh, float *pgh) {
  int ii,mm,nn;
  float g,hh;

  ii = 0;

  for ( nn = 1; nn <= nmax_of_gh; ++nn) {
    for (mm = 0; mm <= nn; ++mm) {
      if (iflag == 1) {
	g = pdata->x;
	hh = pdata->y;
      }
      else {
	g = pdata->z;
	hh = pdata->t;
      }

      pdata++;

      ii++;
      pgh[ii] = g;

      if (mm != 0) {
	ii++;
	pgh[ii] = hh;
      }
    }
  }
}


/****************************************************************************/
/*                                                                          */
/*                           Subroutine extrapsh                            */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Extrapolates linearly a spherical harmonic model with a              */
/*     rate-of-change model.                                                */
/*                                                                          */
/*     Input:                                                               */
/*           date     - date of resulting model (in decimal year)           */
/*           dte1     - date of base model                                  */
/*           nmax1    - maximum degree and order of base model              */
/*           gh1      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of base model                 */
/*           nmax2    - maximum degree and order of rate-of-change model    */
/*           gh2      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of rate-of-change model       */
/*                                                                          */
/*     Output:                                                              */
/*           gha or b - Schmidt quasi-normal internal spherical             */
/*                    harmonic coefficients                                 */
/*           nmax   - maximum degree and order of resulting model           */
/*                                                                          */
/*     FORTRAN                                                              */
/*           A. Zunde                                                       */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 16, 1988                                                */
/*                                                                          */
/****************************************************************************/

static int
extrapsh(struct GeoContext *pctx, float date, float dte1, int nmax1,
	 int nmax2, float *pgh) {
  int   nmax;
  int   k, l;
  int   ii;
  float factor;

  factor = date - dte1;
  if (nmax1 == nmax2) {
    k =  nmax1 * (nmax1 + 2);
    nmax = nmax1;
  }
  else if (nmax1 > nmax2) {
    k = nmax2 * (nmax2 + 2);
    l = nmax1 * (nmax1 + 2);

    for ( ii = k + 1; ii <= l; ++ii)
      pgh[ii] = pctx->gh1[ii];

    nmax = nmax1;
  }
  else {
    k = nmax1 * (nmax1 + 2);
    l = nmax2 * (nmax2 + 2);

    for ( ii = k + 1; ii <= l; ++ii)
      pgh[ii] = factor * pctx->gh2[ii];

    nmax = nmax2;
  }

  for ( ii = 1; ii <= k; ++ii)
    pgh[ii] = pctx->gh1[ii] + factor * pctx->gh2[ii];

   return nmax;
}

/****************************************************************************/
/*                                                                          */
/*                           Subroutine interpsh                            */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Interpolates linearly, in time, between two spherical harmonic       */
/*     models.                                                              */
/*                                                                          */
/*     Input:                                                               */
/*           date     - date of resulting model (in decimal year)           */
/*           dte1     - date of earlier model                               */
/*           nmax1    - maximum degree and order of earlier model           */
/*           gh1      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of earlier model              */
/*           dte2     - date of later model                                 */
/*           nmax2    - maximum degree and order of later model             */
/*           gh2      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of internal model             */
/*                                                                          */
/*     Output:                                                              */
/*           gha or b - coefficients of resulting model                     */
/*           nmax     - maximum degree and order of resulting model         */
/*                                                                          */
/*     FORTRAN                                                              */
/*           A. Zunde                                                       */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 17, 1988                                                */
/*                                                                          */
/****************************************************************************/


static int
interpsh(struct GeoContext *pctx, float date, float dte1, int nmax1,
	 float dte2, int nmax2, float *pgh) {
  int   nmax;
  int   k, l;
  int   ii;
  float factor;

  factor = (date - dte1) / (dte2 - dte1);
  if (nmax1 == nmax2) {
    k =  nmax1 * (nmax1 + 2);
    nmax = nmax1;
  }
  else if (nmax1 > nmax2) {
    k = nmax2 * (nmax2 + 2);
    l = nmax1 * (nmax1 + 2);

    for ( ii = k + 1; ii <= l; ++ii)
      pgh[ii] = pctx->gh1[ii] + factor * (-pctx->gh1[ii]);

    nmax = nmax1;
  }
  else {
    k = nmax1 * (nmax1 + 2);
    l = nmax2 * (nmax2 + 2);
    for ( ii = k + 1; ii <= l; ++ii)
      pgh[ii] = factor * pctx->gh2[ii];

    nmax = nmax2;
  }

  for ( ii = 1; ii <= k; ++ii)
    pgh[ii] = pctx->gh1[ii] + factor * (pctx->gh2[ii] - pctx->gh1[ii]);

   return nmax;
}





/****************************************************************************/
/*                                                                          */
/*                           Subroutine shval3                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Calculates field components from spherical harmonic (sh)             */
/*     models.                                                              */
/*                                                                          */
/*     Input:                                                               */
/*           latitude  - north latitude, in radians                         */
/*           longitude - east longitude, in radians                         */
/*           elev      - elevation above mean sea level                     */
/*           a2,b2     - squares of semi-major and semi-minor axes of       */
/*                       the reference spheroid used for transforming       */
/*                       between geodetic and geocentric coordinates        */
/*                       or components                                      */
/*           nmax      - maximum degree and order of coefficients           */
/*                                                                          */
/*     Output:                                                              */
/*           x         - northward component                                */
/*           y         - eastward component                                 */
/*           z         - vertically-downward component                      */
/*                                                                          */
/*     based on subroutine 'igrf' by D. R. Barraclough and S. R. C. Malin,  */
/*     report no. 71/1, institute of geological sciences, U.K.              */
/*                                                                          */
/*     FORTRAN                                                              */
/*           Norman W. Peddie                                               */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 17, 1988                                                */
/*                                                                          */
/****************************************************************************/


static int
shval3(struct GeoContext *pctx, float flat, float flon, float elev, int nmax, int gh) {
  float earths_radius = 6371.2f;
  float slat;
  float clat;
  float ratio;
  float aa, bb, cc, dd;
  float sd;
  float cd;
  float r;
  float a2;
  float b2;
  float rr = 0.0;
  float fm,fn = 0.0;
  float sl[14];
  float cl[14];
  float p[119];
  float q[119];
  int ii,j,k,l,m,n;
  int npq;
  int ios;
  a2 = 40680925.0f;
  b2 = 40408588.0f;
  ios = 0;
  r = elev;
  slat = (float) sin( flat );
  if (PI / 2 - flat < 0.000017f)
    aa = (float) (PI / 2) - 0.000017f;  /* ~300 ft. from North pole  */
  else if (PI / 2 + flat < 0.000017f)
    aa = (float) -(PI / 2 - 0.000017f); /* ~300 ft. from South pole */
  else
    aa = flat;

  clat = (float) cos( aa );
  sl[1] = (float) sin( flon );
  cl[1] = (float) cos( flon );

  if (gh == 3) {
    pctx->x = 0;
    pctx->y = 0;
    pctx->z = 0;
  }
  else {
    assert(gh == 4);
    pctx->xtemp = 0;
    pctx->ytemp = 0;
    pctx->ztemp = 0;
  }

  sd = 0.0;
  cd = 1.0;
  l = 1;
  n = 0;
  m = 1;
  npq = (nmax * (nmax + 3)) / 2;

  aa = a2 * clat * clat;
  bb = b2 * slat * slat;
  cc = aa + bb;
  dd = (float) sqrt( cc );
  r = (float) sqrt( elev * (elev + 2.0 * dd) + (a2 * aa + b2 * bb) / cc );
  cd = (elev + dd) / r;
  sd = (a2 - b2) / dd * slat * clat / r;
  aa = slat;
  slat = slat * cd - clat * sd;
  clat = clat * cd + aa * sd;

  ratio = earths_radius / r;
  aa = (float) sqrt( 3.0 );
  p[1] = 2.0f * slat;
  p[2] = 2.0f * clat;
  p[3] = 4.5f * slat * slat - 1.5f;
  p[4] = 3.0f * aa * clat * slat;
  q[1] = -clat;
  q[2] = slat;
  q[3] = -3.0f * clat * slat;
  q[4] = aa * (slat * slat - clat * clat);

  for ( k = 1; k <= npq; ++k) {
    if (n < m) {
      m = 0;
      n = n + 1;
      rr = (float) pow(ratio, n + 2);
      fn = (float) n;
    }

    fm = (float) m;
    if (k >= 5) {
      if (m == n) {
	aa = (float) sqrt( 1.0 - 0.5/fm );
	j = k - n - 1;
	p[k] = (1.0f + 1.0f/fm) * aa * clat * p[j];
	q[k] = aa * (clat * q[j] + slat/fm * p[j]);
	sl[m] = sl[m-1] * cl[1] + cl[m-1] * sl[1];
	cl[m] = cl[m-1] * cl[1] - sl[m-1] * sl[1];
      }
      else {
	aa = (float) sqrt( fn*fn - fm*fm );
	bb = (float) sqrt (((fn - 1.0)*(fn-1.0)) - (fm * fm)) / aa;
	cc = (2.0f * fn - 1.0f)/aa;
	ii = k - n;
	j = k - 2 * n + 1;
	p[k] = (fn + 1.0f) * (cc * slat/fn * p[ii] - bb/(fn - 1.0f) * p[j]);
	q[k] = cc * (slat * q[ii] - clat/fn * p[ii]) - bb * q[j];
      }
    }

    if (gh == 3) {
      aa = rr * pctx->gha[l];
    }
    else {
      aa = rr * pctx->ghb[l];
    }
    
    if (m == 0) {
      if (gh == 3) {
	pctx->x += aa * q[k];
	pctx->z -= aa * p[k];
      }
      else {
	pctx->xtemp += aa * q[k];
	pctx->ztemp -= aa * p[k];
      }

      l = l + 1;
    }
    else if (gh == 3) {
      bb = rr * pctx->gha[l+1];
      cc = aa * cl[m] + bb * sl[m];
      pctx->x += cc * q[k];
      pctx->z -= cc * p[k];
      if (clat > 0) {
	pctx->y += (aa * sl[m] - bb * cl[m]) *
	  fm * p[k]/((fn + 1.0f) * clat);
      }
      else {
	pctx->y += (aa * sl[m] - bb * cl[m]) * q[k] * slat;
      }
      l = l + 2;
    }
    else {
      bb = rr * pctx->ghb[l+1];
      cc = aa * cl[m] + bb * sl[m];
      pctx->xtemp += cc * q[k];
      pctx->ztemp -= cc * p[k];
      if (clat > 0) {
	pctx->ytemp += (aa * sl[m] - bb * cl[m]) *
	  fm * p[k]/((fn + 1.0f) * clat);
      }
      else {
	pctx->ytemp += (aa * sl[m] - bb * cl[m]) *
	  q[k] * slat;
      }
      l = l + 2;
    }

    m = m + 1;
  }

  if (gh == 3) {
    aa = pctx->x;
    pctx->x = pctx->x * cd + pctx->z * sd;
    pctx->z = pctx->z * cd - aa * sd;
  }
  else {
    aa = pctx->xtemp;
    pctx->xtemp = pctx->xtemp * cd + pctx->ztemp * sd;
    pctx->ztemp = pctx->ztemp * cd - aa * sd;
  }
  
  return(ios);
}


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

static void
dihf3(struct GeoContext *pctx) {
  int j;
  float sn;
  float h2;
  float hpx;

  float x = pctx->x, y = pctx->y, z = pctx->z;
  float d, i, h, f;

  sn = 0.0001f;

  for (j = 1; j <= 1; ++j) {
    h2 = x*x + y*y;
    h = (float) sqrt(h2);       /* calculate horizontal intensity */
    f = (float) sqrt(h2 + z*z);      /* calculate total intensity */
    if (f < sn) {
      d = NAN;        /* If d and i cannot be determined, */
      i = NAN;        /*       set equal to 999.0         */
    }
    else {
      i = (float) atan2(z, h);
      if (h < sn) {
	d = NAN;
      }
      else {
	hpx = h + x;
	if (hpx < sn)
	  d = PI;
	else
	  d = 2.0f * (float) atan2(y, hpx);
      }
    }
  }

  pctx->d = d;
  pctx->i = i;
  pctx->h = h;
  pctx->f = f;
}

static void
dihf4(struct GeoContext *pctx) {
  int j;
  float sn;
  float h2;
  float hpx;

  float xtemp = pctx->xtemp, ytemp = pctx->ytemp, ztemp = pctx->ztemp;
  float dtemp, itemp, htemp, ftemp;

  sn = 0.0001f;

  for (j = 1; j <= 1; ++j) {
    h2 = xtemp*xtemp + ytemp*ytemp;
    htemp = (float) sqrt(h2);
    ftemp = (float) sqrt(h2 + ztemp*ztemp);
    if (ftemp < sn) {
      dtemp = NAN;    /* If d and i cannot be determined, */
      itemp = NAN;    /*       set equal to 999.0         */
    }
    else {
      itemp = (float) atan2(ztemp, htemp);
      if (htemp < sn) {
	dtemp = NAN;
      }
      else {
	hpx = htemp + xtemp;
	if (hpx < sn)
	  dtemp = PI;
	else
	  dtemp = 2.0f * (float) atan2(ytemp, hpx);
      }
    }
  }

  pctx->dtemp = dtemp;
  pctx->itemp = itemp;
  pctx->htemp = htemp;
  pctx->ftemp = ftemp;
}
