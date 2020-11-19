/*  worldpos.c -- WCS Algorithms from Classic AIPS.
 *  June 4, 1997
 *  Copyright (C) 1994
 *  Associated Universities, Inc. Washington DC, USA.
 *  With code added by Doug Mink, Smithsonian Astrophysical Observatory

 * Module:	worldpos.c
 * Purpose:	Perform forward and reverse WCS computations for 8 projections
 * Subroutine:	worldpos() converts from pixel location to RA,Dec 
 * Subroutine:	worldpix() converts from RA,Dec         to pixel location   

    This library is free software; you can redistribute it and/or modify it
    under the terms of the GNU Library General Public License as published by
    the Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.
   
    This library is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
    License for more details.
   
    You should have received a copy of the GNU Library General Public License
    along with this library; if not, write to the Free Software Foundation,
    Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
   
    Correspondence concerning AIPS should be addressed as follows:
	   Internet email: aipsmail@nrao.edu
	   Postal address: AIPS Group
	                   National Radio Astronomy Observatory
	                   520 Edgemont Road
	                   Charlottesville, VA 22903-2475 USA

	         -=-=-=-=-=-=-

    These two ANSI C functions, worldpos() and worldpix(), perform
    forward and reverse WCS computations for 8 types of projective
    geometries ("-SIN", "-TAN", "-ARC", "-NCP", "-GLS", "-MER", "-AIT"
    and "-STG"):

	worldpos() converts from pixel location to RA,Dec 
	worldpix() converts from RA,Dec         to pixel location   

    where "(RA,Dec)" are more generically (long,lat). These functions
    are based on the WCS implementation of Classic AIPS, an
    implementation which has been in production use for more than ten
    years. See the two memos by Eric Greisen

	ftp://fits.cv.nrao.edu/fits/documents/wcs/aips27.ps.Z
	ftp://fits.cv.nrao.edu/fits/documents/wcs/aips46.ps.Z

    for descriptions of the 8 projective geometries and the
    algorithms.  Footnotes in these two documents describe the
    differences between these algorithms and the 1993-94 WCS draft
    proposal (see URL below). In particular, these algorithms support
    ordinary field rotation, but not skew geometries (CD or PC matrix
    cases). Also, the MER and AIT algorithms work correctly only for
    CRVALi=(0,0). Users should note that GLS projections with yref!=0
    will behave differently in this code than in the draft WCS
    proposal.  The NCP projection is now obsolete (it is a special
    case of SIN).  WCS syntax and semantics for various advanced
    features is discussed in the draft WCS proposal by Greisen and
    Calabretta at:
    
	ftp://fits.cv.nrao.edu/fits/documents/wcs/wcs.all.ps.Z
    
	        -=-=-=-

    The original version of this code was Emailed to D.Wells on
    Friday, 23 September by Bill Cotton <bcotton@gorilla.cv.nrao.edu>,
    who described it as a "..more or less.. exact translation from the
    AIPSish..". Changes were made by Don Wells <dwells@nrao.edu>
    during the period October 11-13, 1994:
    1) added GNU license and header comments
    2) added testpos.c program to perform extensive circularity tests
    3) changed float-->double to get more than 7 significant figures
    4) testpos.c circularity test failed on MER and AIT. B.Cotton
       found that "..there were a couple of lines of code [in] the wrong
       place as a result of merging several Fortran routines." 
    5) testpos.c found 0h wraparound in worldpix() and worldpos().
    6) E.Greisen recommended removal of various redundant if-statements,
       and addition of a 360d difference test to MER case of worldpos(). 
    7) D.Mink changed input to data structure and implemented rotation matrix.
*/
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "wcs.h"

int
worldpos (xpix, ypix, wcs, xpos, ypos)

/* Routine to determine accurate position for pixel coordinates */
/* returns 0 if successful otherwise 1 = angle too large for projection; */
/* does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections */
/* anything else is linear */

/* Input: */
double	xpix;		/* x pixel number  (RA or long without rotation) */
double	ypix;		/* y pixel number  (dec or lat without rotation) */
struct WorldCoor *wcs;		/* WCS parameter structure */

/* Output: */
double	*xpos;		/* x (RA) coordinate (deg) */
double	*ypos;		/* y (dec) coordinate (deg) */

{
  double cosr, sinr, dx, dy, dz, temp;
  double sins, coss, dect, rat, dt, l, m, mg, da, dd, cos0, sin0;
  double dec0, ra0, decout, raout;
  double geo1, geo2, geo3;
  double cond2r=1.745329252e-2;
  double twopi = 6.28318530717959;
  double deps = 1.0e-5;

  /* Structure elements */
  double xref;		/* x reference coordinate value (deg) */
  double yref;		/* y reference coordinate value (deg) */
  double xrefpix;	/* x reference pixel */
  double yrefpix;	/* y reference pixel */
  double xinc;		/* x coordinate increment (deg) */
  double yinc;		/* y coordinate increment (deg) */
  double rot;		/* rotation (deg)  (from N through E) */
  int itype = wcs->pcode;

/* Set local projection parameters */
  xref = wcs->xref;
  yref = wcs->yref;
  xrefpix = wcs->xrefpix;
  yrefpix = wcs->yrefpix;
  xinc = wcs->xinc;
  yinc = wcs->yinc;

/* Offset from ref pixel */
  dx = xpix - xrefpix;
  dy = ypix - yrefpix;

/* Scale and rotate using CD matrix */
  if (wcs->rotmat) {
    temp = dx * wcs->cd11 + dy * wcs->cd12;
    dy = dx * wcs->cd21 + dy * wcs->cd22;
    dx = temp;
    }
  else {

/* Check axis increments - bail out if either 0 */
    if ((xinc==0.0) || (yinc==0.0)) {
      *xpos=0.0;
      *ypos=0.0;
      return 2;
      }

/* Scale using CDELT */
    dx = dx * xinc;
    dy = dy * yinc;

/* Take out rotation from CROTA */
    rot = wcs->rot;
    cosr = wcs->crot;
    sinr = wcs->srot;
    if (rot != 0.0) {
      temp = dx * cosr - dy * sinr;
      dy = dy * cosr + dx * sinr;
      dx = temp;
      }
    }

/* Default, linear result for error or pixel return  */
  *xpos = xref + dx;
  *ypos = yref + dy;
  if (itype < 0)
    return 0;

/* Convert to radians  */
  if (wcs->coorflip) {
    dec0 = degrad (xref);
    ra0 = degrad (yref);
    temp = dx;
    dx = dy;
    dy = temp;
    }
  else {
    ra0 = degrad (xref);
    dec0 = degrad (yref);
    }
  l = degrad (dx);
  m = degrad (dy);
  sins = l*l + m*m;
  decout = 0.0;
  raout = 0.0;
  cos0 = cos(dec0);
  sin0 = sin(dec0);

/* process by case  */
  switch (itype) {
    case -1:   /* pixel */
    case 0:   /* linear */
      rat =  ra0 + l;
      dect = dec0 + m;
      break;
    case 1:   /* -SIN sin*/ 
      if (sins>1.0) return 1;
      coss = sqrt (1.0 - sins);
      dt = sin0 * coss + cos0 * m;
      if ((dt>1.0) || (dt<-1.0)) return 1;
      dect = asin (dt);
      rat = cos0 * coss - sin0 * m;
      if ((rat==0.0) && (l==0.0)) return 1;
      rat = atan2 (l, rat) + ra0;
      break;
    case 2:   /* -TAN tan */
      if (sins>1.0) return 1;
      dect = cos0 - m * sin0;
      if (dect==0.0) return 1;
      rat = ra0 + atan2 (l, dect);
      dect = atan (cos(rat-ra0) * (m * cos0 + sin0) / dect);
      break;
    case 3:   /* -ARC Arc*/
      if (sins>=twopi*twopi/4.0) return 1;
      sins = sqrt(sins);
      coss = cos (sins);
      if (sins!=0.0) sins = sin (sins) / sins;
      else
	sins = 1.0;
      dt = m * cos0 * sins + sin0 * coss;
      if ((dt>1.0) || (dt<-1.0)) return 1;
      dect = asin (dt);
      da = coss - dt * sin0;
      dt = l * sins * cos0;
      if ((da==0.0) && (dt==0.0)) return 1;
      rat = ra0 + atan2 (dt, da);
      break;
    case 4:   /* -NCP North celestial pole*/
      dect = cos0 - m * sin0;
      if (dect==0.0) return 1;
      rat = ra0 + atan2 (l, dect);
      dt = cos (rat-ra0);
      if (dt==0.0) return 1;
      dect = dect / dt;
      if ((dect>1.0) || (dect<-1.0)) return 1;
      dect = acos (dect);
      if (dec0<0.0) dect = -dect;
      break;
    case 5:   /* -GLS global sinusoid */
      dect = dec0 + m;
      if (fabs(dect)>twopi/4.0) return 1;
      coss = cos (dect);
      if (fabs(l)>twopi*coss/2.0) return 1;
      rat = ra0;
      if (coss>deps) rat = rat + l / coss;
      break;
    case 6:   /* -MER mercator*/
      dt = yinc * cosr + xinc * sinr;
      if (dt==0.0) dt = 1.0;
      dy = degrad (yref/2.0 + 45.0);
      dx = dy + dt / 2.0 * cond2r;
      dy = log (tan (dy));
      dx = log (tan (dx));
      geo2 = degrad (dt) / (dx - dy);
      geo3 = geo2 * dy;
      geo1 = cos (degrad (yref));
      if (geo1<=0.0) geo1 = 1.0;
      rat = l / geo1 + ra0;
      if (fabs(rat - ra0) > twopi) return 1; /* added 10/13/94 DCW/EWG */
      dt = 0.0;
      if (geo2!=0.0) dt = (m + geo3) / geo2;
      dt = exp (dt);
      dect = 2.0 * atan (dt) - twopi / 4.0;
      break;
    case 7:   /* -AIT Aitoff*/
      dt = yinc*cosr + xinc*sinr;
      if (dt==0.0) dt = 1.0;
      dt = degrad (dt);
      dy = degrad (yref);
      dx = sin(dy+dt)/sqrt((1.0+cos(dy+dt))/2.0) -
	  sin(dy)/sqrt((1.0+cos(dy))/2.0);
      if (dx==0.0) dx = 1.0;
      geo2 = dt / dx;
      dt = xinc*cosr - yinc* sinr;
      if (dt==0.0) dt = 1.0;
      dt = degrad (dt);
      dx = 2.0 * cos(dy) * sin(dt/2.0);
      if (dx==0.0) dx = 1.0;
      geo1 = dt * sqrt((1.0+cos(dy)*cos(dt/2.0))/2.0) / dx;
      geo3 = geo2 * sin(dy) / sqrt((1.0+cos(dy))/2.0);
      rat = ra0;
      dect = dec0;
      if ((l==0.0) && (m==0.0)) break;
      dz = 4.0 - l*l/(4.0*geo1*geo1) - ((m+geo3)/geo2)*((m+geo3)/geo2) ;
      if ((dz>4.0) || (dz<2.0)) return 1;;
      dz = 0.5 * sqrt (dz);
      dd = (m+geo3) * dz / geo2;
      if (fabs(dd)>1.0) return 1;;
      dd = asin (dd);
      if (fabs(cos(dd))<deps) return 1;;
      da = l * dz / (2.0 * geo1 * cos(dd));
      if (fabs(da)>1.0) return 1;;
      da = asin (da);
      rat = ra0 + 2.0 * da;
      dect = dd;
      break;
    case 8:   /* -STG Sterographic*/
      dz = (4.0 - sins) / (4.0 + sins);
      if (fabs(dz)>1.0) return 1;
      dect = dz * sin0 + m * cos0 * (1.0+dz) / 2.0;
      if (fabs(dect)>1.0) return 1;
      dect = asin (dect);
      rat = cos(dect);
      if (fabs(rat)<deps) return 1;
      rat = l * (1.0+dz) / (2.0 * rat);
      if (fabs(rat)>1.0) return 1;
      rat = asin (rat);
      mg = 1.0 + sin(dect) * sin0 + cos(dect) * cos0 * cos(rat);
      if (fabs(mg)<deps) return 1;
      mg = 2.0 * (sin(dect) * cos0 - cos(dect) * sin0 * cos(rat)) / mg;
      if (fabs(mg-m)>deps) rat = twopi/2.0 - rat;
      rat = ra0 + rat;
      break;
  }
/*  return ra in range  */
  raout = rat;
  decout = dect;
  if (raout-ra0>twopi/2.0) raout = raout - twopi;
  if (raout-ra0<-twopi/2.0) raout = raout + twopi;
  if (raout < 0.0) raout += twopi; /* added by DCW 10/12/94 */

/*  correct units back to degrees  */
  *xpos = raddeg (raout);
  *ypos = raddeg (decout);

  return 0;
}  /* End of worldpos */


int
worldpix (xpos, ypos, wcs, xpix, ypix)

/*-----------------------------------------------------------------------*/
/* routine to determine accurate pixel coordinates for an RA and Dec     */
/* returns 0 if successful otherwise:                                    */
/*  1 = angle too large for projection;                                  */
/*  2 = bad values                                                       */
/* does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections            */
/* anything else is linear                                               */

/* Input: */
double	xpos;		/* x (RA) coordinate (deg) */
double	ypos;		/* y (dec) coordinate (deg) */
struct WorldCoor *wcs;	/* WCS parameter structure */

/* Output: */
double	*xpix;		/* x pixel number  (RA or long without rotation) */
double	*ypix;		/* y pixel number  (dec or lat without rotation) */
{
  double dx, dy, ra0, dec0, ra, dec, coss, sins, dt, da, dd, sint;
  double l, m, geo1, geo2, geo3, sinr, cosr, temp;
  double cond2r=1.745329252e-2, deps=1.0e-5, twopi=6.28318530717959;

/* Structure elements */
  double xref;		/* x reference coordinate value (deg) */
  double yref;		/* y reference coordinate value (deg) */
  double xrefpix;	/* x reference pixel */
  double yrefpix;	/* y reference pixel */
  double xinc;		/* x coordinate increment (deg) */
  double yinc;		/* y coordinate increment (deg) */
  double rot;		/* rotation (deg)  (from N through E) */
  int itype;

/* Set local projection parameters */
  xref = wcs->xref;
  yref = wcs->yref;
  xrefpix = wcs->xrefpix;
  yrefpix = wcs->yrefpix;
  xinc = wcs->xinc;
  yinc = wcs->yinc;
  rot = wcs->rot;
  cosr = wcs->crot;
  sinr = wcs->srot;

/* Projection type */
  itype = wcs->pcode;

  /* 0h wrap-around tests added by D.Wells 10/12/94: */
  dt = (xpos - xref);
  if (itype >= 0) {
    if (dt > 180.0) xpos -= 360.0;
    if (dt < -180.0) xpos += 360.0;
    /* NOTE: changing input argument xpos is OK (call-by-value in C!) */
    }

/* Nonlinear position */
  if (itype > 0 && itype < 9) {
    if (wcs->coorflip) {
      dec0 = degrad (xref);
      ra0 = degrad (yref);
      }
    else {
      ra0 = degrad (xref);
      dec0 = degrad (yref);
      }
    ra = degrad (xpos);
    dec = degrad (ypos);

/* Compute direction cosine */
    coss = cos (dec);
    sins = sin (dec);
    l = sin(ra-ra0) * coss;
    sint = sins * sin(dec0) + coss * cos(dec0) * cos(ra-ra0);
    }

/* Process by case  */
  switch (itype) {
    case 1:   /* -SIN sin*/ 
	 if (sint<0.0) return 1;
	 m = sins * cos(dec0) - coss * sin(dec0) * cos(ra-ra0);
      break;
    case 2:   /* -TAN tan */
	 if (sint<=0.0) return 1;
 	 m = sins * sin(dec0) + coss * cos(dec0) * cos(ra-ra0);
	 l = l / m;
	 m = (sins * cos(dec0) - coss * sin(dec0) * cos(ra-ra0)) / m;
      break;
    case 3:   /* -ARC Arc*/
	 m = sins * sin(dec0) + coss * cos(dec0) * cos(ra-ra0);
	 if (m<-1.0) m = -1.0;
	 if (m>1.0) m = 1.0;
	 m = acos (m);
	 if (m!=0) 
	    m = m / sin(m);
	 else
	    m = 1.0;
	 l = l * m;
	 m = (sins * cos(dec0) - coss * sin(dec0) * cos(ra-ra0)) * m;
      break;
    case 4:   /* -NCP North celestial pole*/
	 if (dec0==0.0) 
	     return 1;  /* can't stand the equator */
	 else
	   m = (cos(dec0) - coss * cos(ra-ra0)) / sin(dec0);
      break;
    case 5:   /* -GLS global sinusoid */
	 dt = ra - ra0;
	 if (fabs(dec)>twopi/4.0) return 1;
	 if (fabs(dec0)>twopi/4.0) return 1;
	 m = dec - dec0;
	 l = dt * coss;
      break;
    case 6:   /* -MER mercator*/
	 dt = yinc * cosr + xinc * sinr;
	 if (dt==0.0) dt = 1.0;
	 dy = degrad (yref/2.0 + 45.0);
	 dx = dy + dt / 2.0 * cond2r;
	 dy = log (tan (dy));
	 dx = log (tan (dx));
	 geo2 = degrad (dt) / (dx - dy);
	 geo3 = geo2 * dy;
	 geo1 = cos (degrad (yref));
	 if (geo1<=0.0) geo1 = 1.0;
	 dt = ra - ra0;
	 l = geo1 * dt;
	 dt = dec / 2.0 + twopi / 8.0;
	 dt = tan (dt);
	 if (dt<deps) return 2;
	 m = geo2 * log (dt) - geo3;
	 break;
    case 7:   /* -AIT Aitoff*/
	 l = 0.0;
	 m = 0.0;
	 da = (ra - ra0) / 2.0;
	 if (fabs(da)>twopi/4.0) return 1;
	 dt = yinc*cosr + xinc*sinr;
	 if (dt==0.0) dt = 1.0;
	 dt = degrad (dt);
	 dy = degrad (yref);
	 dx = sin(dy+dt)/sqrt((1.0+cos(dy+dt))/2.0) -
	     sin(dy)/sqrt((1.0+cos(dy))/2.0);
	 if (dx==0.0) dx = 1.0;
	 geo2 = dt / dx;
	 dt = xinc*cosr - yinc* sinr;
	 if (dt==0.0) dt = 1.0;
	 dt = degrad (dt);
	 dx = 2.0 * cos(dy) * sin(dt/2.0);
	 if (dx==0.0) dx = 1.0;
	 geo1 = dt * sqrt((1.0+cos(dy)*cos(dt/2.0))/2.0) / dx;
	 geo3 = geo2 * sin(dy) / sqrt((1.0+cos(dy))/2.0);
	 dt = sqrt ((1.0 + cos(dec) * cos(da))/2.0);
	 if (fabs(dt)<deps) return 3;
	 l = 2.0 * geo1 * cos(dec) * sin(da) / dt;
	 m = geo2 * sin(dec) / dt - geo3;
      break;
    case 8:   /* -STG Sterographic*/
	 da = ra - ra0;
	 if (fabs(dec)>twopi/4.0) return 1;
	 dd = 1.0 + sins * sin(dec0) + coss * cos(dec0) * cos(da);
	 if (fabs(dd)<deps) return 1;
	 dd = 2.0 / dd;
	 l = l * dd;
	 m = dd * (sins * cos(dec0) - coss * sin(dec0) * cos(da));
      break;
  }  /* end of itype switch */

/* Back to degrees  */
  if (itype > 0 && itype < 9) {
    dx = raddeg (l);
    dy = raddeg (m);
    }

/* For linear or pixel projection */
  else {
    dx = xpos - xref;
    dy = ypos - yref;
    }

/* Scale and rotate using CD matrix */
  if (wcs->rotmat) {
    temp = dx * wcs->dc11 + dy * wcs->dc12;
    dy = dx * wcs->dc21 + dy * wcs->dc22;
    dx = temp;
    }
  else {

/* Correct for rotation */
    if (rot!=0.0) {
      temp = dx*cosr + dy*sinr;
      dy = dy*cosr - dx*sinr;
      dx = temp;
      }

/* Scale using CDELT */
    if (xinc != 0.)
      dx = dx / xinc;
    if (yinc != 0.)
      dy = dy / yinc;
    }

  if (wcs->coorflip) {
    temp = dx;
    dx = dy;
    dy = temp;
    }

/* Convert to pixels  */
  *xpix = dx + xrefpix;
  *ypix = dy + yrefpix;

  return 0;
}  /* end worldpix */

/* Oct 26 1995	Fix bug which interchanged RA and Dec twice when coorflip
 * Oct 31 1996	Fix CD matrix use in WORLDPIX
 * Nov  4 1996	Eliminate extra code for linear projection in WORLDPIX
 * Nov  5 1996	Add coordinate flip in WORLDPIX
 * May 22 1997	Avoid angle wraparound when CTYPE is pixel
 * Jun  4 1997	Return without angle conversion from worldpos if type is PIXEL
 */
