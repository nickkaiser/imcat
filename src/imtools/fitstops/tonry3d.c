/*  -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "myf2c.h"
/* #include "stdlib.h" */

/* Common Block Declarations */

struct {
    real a1, a2, a3, b1, b2, b3, b4;
} mgoplt3b_;

#define mgoplt3b_1 mgoplt3b_

struct {
    real xx[2000], yy[2000];
    integer kk, ll;
} mgonxtv1_;

#define mgonxtv1_1 mgonxtv1_

/* Table of constant values */

static integer c__1 = 1;

/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*/
/*C                  Mongo Interactive Graphics Software                    CC
*/
/*C                 Copyright (c) 1987, 1994 - John Tonry.                  CC
*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*/
/* *** argument list *** */
/* 	SUBROUTINE MGOPLT3D(A,M,N,WORK,ALT,AZ,ZFAC,ZOFF) */

/* 	A	- REAL DATA ARRAY, REPRESENTS HEIGHT OF SURFACE AS */
/* 		  FUNCTION MGOOF LOCATION IN PLANE */
/* 	M,N	- DIMENSIONS OF DATA ARRAY A */
/* 	WORK	- REAL WORK ARRAY, LENGTH .GE. 4*MIN(M,N) */
/* 	ALT,AZ	- ALTITUDE,AZIMUTH VIEWING ANGLES IN DEGREES */
/* 	ZFAC	- SCALING OF Z-AXIS (INVERSE OF RANGE, I.E. DEFAULT 1/DATAMAX) */
/* 	ZOFF	- OFFSET OF Z-ORIGIN IN DATA UNITS */

/* Note: A(I,J) altered to A(M-I+1,J) to correct x parity flip */

/* Subroutine */ 
int mgoplt3d_(float *a, int *m, int *n, float *alt, float *az, float *zfac, float *zoff, float width, float height)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Builtin functions */
    double sin(), cos();

    /* Local variables */
    static integer ibeg, jbeg, lmax;
    static real xoff, yoff, talt, xlen, ylen;
    static integer lnth;
    extern /* Subroutine */ int mgonxtvu_();
    static integer i__, j, ilast, jlast, istep, jstep, ic, ll;
    static real xx, yy;
    static integer ifirst, jfirst;
    static real gx1, gx2, gy1, gy2, cal, caz, sal;
    static integer lli, ier, iaz;
    static real saz, taz, xsc, ysc, zfactor;

	float *work;
	work = (float *) calloc(4 * *n, sizeof(float)); 

/*      INCLUDE 'mongo.par' */
    /* Parameter adjustments */
    a_dim1 = *m;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --work;

    /* Function Body */
    gx1 = (float)0.;
    gx2 = width;
    gy1 = (float)0.;
    gy2 = height;

    xlen = (gx2 - gx1) * 2 / 3;
    ylen = xlen;
    xoff = (gx2 + gx1) * (float).5;
    yoff = (gy2 + gy1) * (float).5;
    zfactor = *zfac * (gy2 - gy1) / 4;
    lmax = min(*m,*n) << 1;
    taz = *az * (float).0174532925;
    talt = *alt * (float).0174532925;
    saz = sin(taz);
    caz = cos(taz);
    sal = sin(talt);
    cal = cos(talt);
    xsc = xlen / (real) (*n - 1);
    ysc = ylen / (real) (*m - 1);
    mgoplt3b_1.a1 = caz * xsc;
    mgoplt3b_1.a2 = -saz * ysc;
    mgoplt3b_1.a3 = xoff - (mgoplt3b_1.a1 * (real) (*n + 1) + mgoplt3b_1.a2 * 
	    (real) (*m + 1)) * (float).5;
    mgoplt3b_1.b1 = saz * sal * xsc;
    mgoplt3b_1.b2 = caz * sal * ysc;
    mgoplt3b_1.b3 = zfactor * cal;
    mgoplt3b_1.b4 = mgoplt3b_1.b3 * *zoff + yoff - (mgoplt3b_1.b1 * (real) (*
	    n + 1) + mgoplt3b_1.b2 * (real) (*m + 1)) * (float).5;
    iaz = 1;
    if (mgoplt3b_1.a1 <= (float)0.) {
	++iaz;
    }
    if (mgoplt3b_1.a2 <= (float)0.) {
	iaz += 2;
    }
    switch ((int)iaz) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L10;
	case 4:  goto L20;
    }
L10:
    ifirst = 1;
    istep = 1;
    ilast = *m;
    goto L30;
L20:
    ifirst = *m;
    istep = -1;
    ilast = 1;
L30:
    switch ((int)iaz) {
	case 1:  goto L50;
	case 2:  goto L50;
	case 3:  goto L40;
	case 4:  goto L40;
    }
L40:
    jfirst = 1;
    jstep = 1;
    jlast = *n;
    goto L60;
L50:
    jfirst = *n;
    jstep = -1;
    jlast = 1;
L60:
    switch ((int)iaz) {
	case 1:  goto L64;
	case 2:  goto L62;
	case 3:  goto L62;
	case 4:  goto L64;
    }
L62:
    lli = 1;
    goto L66;
L64:
    lli = -1;
L66:
    ic = 0;
    ibeg = ifirst + istep;
L70:
/* Computing MIN */
    i__2 = ((i__1 = ibeg - ifirst, abs(i__1)) << 1) + 1;
    lnth = min(i__2,lmax);
    if (lli == -1) {
	goto L72;
    }
    ll = 0;
    goto L74;
L72:
    ll = lnth + 1;
L74:
    i__ = ibeg;
    j = jfirst;
    xx = (real) j;
    yy = (real) i__;
    ll += lli;
    work[ll] = mgoplt3b_1.a1 * xx + mgoplt3b_1.a2 * yy + mgoplt3b_1.a3;
    work[ll + lmax] = mgoplt3b_1.b1 * xx + mgoplt3b_1.b2 * yy + mgoplt3b_1.b3 
	    * (a[*m - i__ + 1 + j * a_dim1] + *zoff) + mgoplt3b_1.b4;
L80:
    i__ -= istep;
    yy = (real) i__;
    ll += lli;
    work[ll] = mgoplt3b_1.a1 * xx + mgoplt3b_1.a2 * yy + mgoplt3b_1.a3;
    work[ll + lmax] = mgoplt3b_1.b1 * xx + mgoplt3b_1.b2 * yy + mgoplt3b_1.b3 
	    * (a[*m - i__ + 1 + j * a_dim1] + *zoff) + mgoplt3b_1.b4;
    if (j == jlast) {
	goto L85;
    }
    j += jstep;
    xx = (real) j;
    ll += lli;
    work[ll] = mgoplt3b_1.a1 * xx + mgoplt3b_1.a2 * yy + mgoplt3b_1.a3;
    work[ll + lmax] = mgoplt3b_1.b1 * xx + mgoplt3b_1.b2 * yy + mgoplt3b_1.b3 
	    * (a[*m - i__ + 1 + j * a_dim1] + *zoff) + mgoplt3b_1.b4;
    if (i__ != ifirst) {
	goto L80;
    }
L85:
    mgonxtvu_(&ic, &work[1], &work[lmax + 1], &lnth, &ier);
    if (ier != 0) {
	return 0;
    }
    ic = 1;
    if (ibeg == ilast) {
	goto L90;
    }
    ibeg += istep;
    goto L70;
L90:
    jbeg = jfirst;
L100:
/* Computing MIN */
    i__2 = ((i__1 = jbeg - jlast, abs(i__1)) << 1) + 1;
    lnth = min(i__2,lmax);
    if (lli == -1) {
	goto L102;
    }
    ll = 0;
    goto L104;
L102:
    ll = lnth + 1;
L104:
    i__ = ilast;
    j = jbeg;
    xx = (real) j;
    yy = (real) i__;
    ll += lli;
    work[ll] = mgoplt3b_1.a1 * xx + mgoplt3b_1.a2 * yy + mgoplt3b_1.a3;
    work[ll + lmax] = mgoplt3b_1.b1 * xx + mgoplt3b_1.b2 * yy + mgoplt3b_1.b3 
	    * (a[*m - i__ + 1 + j * a_dim1] + *zoff) + mgoplt3b_1.b4;
L110:
    j += jstep;
    xx = (real) j;
    ll += lli;
    work[ll] = mgoplt3b_1.a1 * xx + mgoplt3b_1.a2 * yy + mgoplt3b_1.a3;
    work[ll + lmax] = mgoplt3b_1.b1 * xx + mgoplt3b_1.b2 * yy + mgoplt3b_1.b3 
	    * (a[*m - i__ + 1 + j * a_dim1] + *zoff) + mgoplt3b_1.b4;
    if (i__ == ifirst) {
	goto L120;
    }
    i__ -= istep;
    yy = (real) i__;
    ll += lli;
    work[ll] = mgoplt3b_1.a1 * xx + mgoplt3b_1.a2 * yy + mgoplt3b_1.a3;
    work[ll + lmax] = mgoplt3b_1.b1 * xx + mgoplt3b_1.b2 * yy + mgoplt3b_1.b3 
	    * (a[*m - i__ + 1 + j * a_dim1] + *zoff) + mgoplt3b_1.b4;
    if (j != jlast) {
	goto L110;
    }
L120:
    mgonxtvu_(&c__1, &work[1], &work[lmax + 1], &lnth, &ier);
    if (ier != 0) {
	return 0;
    }
    jbeg += jstep;
    if (jbeg == jlast) {
	return 0;
    }
    goto L100;
} /* mgoplt3d_ */

/* Subroutine */ int mgonxtvu_(ic, x, y, n, ier)
integer *ic;
real *x, *y;
integer *n, *ier;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real frac;
    static integer i__, j;
    static real x0, x1;
    static integer ii, jj;
    static real xi, yi, xl, yl, px, py, ya0, yb0, ya1, yb1;
    static integer isw;
    extern doublereal mgoalin_();
    extern /* Subroutine */ int mgoline_(), mgooutp_();
    static integer iov0, iov1;

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    if (*ic != 0) {
	goto L20;
    }
    if (*n > 2000) {
	goto L500;
    }
    mgonxtv1_1.ll = 2000 - *n + 1;
    i__ = mgonxtv1_1.ll;
    mgonxtv1_1.xx[i__ - 1] = x[1];
    mgonxtv1_1.yy[i__ - 1] = y[1];
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	++i__;
	mgonxtv1_1.xx[i__ - 1] = x[j];
	mgonxtv1_1.yy[i__ - 1] = y[j];
	mgoline_(&mgonxtv1_1.xx[i__ - 2], &mgonxtv1_1.yy[i__ - 2], &
		mgonxtv1_1.xx[i__ - 1], &mgonxtv1_1.yy[i__ - 1]);
/* L10: */
    }
    *ier = 0;
    return 0;
L20:
    if (*ier != 0) {
	return 0;
    }
    ii = 1;
    jj = mgonxtv1_1.ll;
    mgonxtv1_1.kk = 0;
    ya0 = y[1];
    yb0 = mgonxtv1_1.yy[mgonxtv1_1.ll - 1];
    if (x[1] - mgonxtv1_1.xx[mgonxtv1_1.ll - 1] <= (float)0.) {
	goto L30;
    } else {
	goto L70;
    }
L30:
    px = x[1];
    py = ya0;
L40:
    mgooutp_(&x[ii], &y[ii], ier);
    if (ii == *n) {
	goto L360;
    }
    ++ii;
    ya0 = y[ii];
    if (x[ii] > mgonxtv1_1.xx[mgonxtv1_1.ll - 1]) {
	goto L50;
    }
    mgoline_(&px, &py, &x[ii], &ya0);
    px = x[ii];
    py = ya0;
    goto L40;
L50:
    --ii;
    xl = x[ii];
    yl = y[ii];
    ya0 = mgoalin_(&x[ii], &x[ii + 1], &y[ii], &y[ii + 1], &mgonxtv1_1.xx[
	    mgonxtv1_1.ll - 1]);
    x0 = mgonxtv1_1.xx[mgonxtv1_1.ll - 1];
    if (ya0 > yb0) {
	goto L90;
    }
    mgoline_(&px, &py, &x0, &ya0);
    px = x0;
    py = ya0;
    mgooutp_(&x0, &ya0, ier);
    mgooutp_(&x0, &yb0, ier);
    goto L100;
L70:
    mgooutp_(&mgonxtv1_1.xx[jj - 1], &mgonxtv1_1.yy[jj - 1], ier);
    if (jj == 2000) {
	goto L380;
    }
    ++jj;
    yb0 = mgonxtv1_1.yy[jj - 1];
    if (x[1] - mgonxtv1_1.xx[jj - 1] >= (float)0.) {
	goto L70;
    } else {
	goto L80;
    }
L80:
    --jj;
    yb0 = mgoalin_(&mgonxtv1_1.xx[jj - 1], &mgonxtv1_1.xx[jj], &mgonxtv1_1.yy[
	    jj - 1], &mgonxtv1_1.yy[jj], &x[1]);
    x0 = x[1];
    if (ya0 <= yb0) {
	goto L100;
    }
    mgooutp_(&x0, &yb0, ier);
    mgooutp_(&x0, &ya0, ier);
    xl = x0;
    yl = ya0;
L90:
    iov0 = 1;
    goto L120;
L100:
    iov0 = 0;
L120:
    if (ii == *n) {
	goto L300;
    }
    if (jj == 2000) {
	goto L310;
    }
    if (x[ii + 1] > mgonxtv1_1.xx[jj]) {
	goto L130;
    }
    isw = 1;
    ++ii;
    x1 = x[ii];
    ya1 = y[ii];
    yb1 = mgoalin_(&mgonxtv1_1.xx[jj - 1], &mgonxtv1_1.xx[jj], &mgonxtv1_1.yy[
	    jj - 1], &mgonxtv1_1.yy[jj], &x1);
    goto L140;
L130:
    if (mgonxtv1_1.xx[jj] >= x[*n]) {
	goto L340;
    }
    isw = -1;
    ++jj;
    x1 = mgonxtv1_1.xx[jj - 1];
    ya1 = mgoalin_(&x[ii], &x[ii + 1], &y[ii], &y[ii + 1], &x1);
    yb1 = mgonxtv1_1.yy[jj - 1];
L140:
    if (ya1 <= yb1) {
	goto L160;
    }
    iov1 = 1;
    if (iov0 == 0) {
	goto L170;
    }
L150:
    if (isw == -1) {
	goto L200;
    }
    mgooutp_(&x1, &ya1, ier);
    mgoline_(&xl, &yl, &x1, &ya1);
    xl = x1;
    yl = ya1;
    goto L200;
L160:
    iov1 = 0;
    if (iov0 == 0) {
	goto L190;
    }
L170:
    frac = (yb0 - ya0) / (ya1 - yb1 + yb0 - ya0);
    xi = (x1 - x0) * frac + x0;
    yi = (ya1 - ya0) * frac + ya0;
    mgooutp_(&xi, &yi, ier);
    if (iov0 == 0) {
	goto L180;
    }
    mgoline_(&xl, &yl, &xi, &yi);
    xl = xi;
    yl = yi;
    goto L190;
L180:
    xl = xi;
    yl = yi;
    goto L150;
L190:
    if (isw == 1) {
	goto L200;
    }
    mgooutp_(&mgonxtv1_1.xx[jj - 1], &mgonxtv1_1.yy[jj - 1], ier);
L200:
    if (*ier != 0) {
	return 0;
    }
    x0 = x1;
    ya0 = ya1;
    yb0 = yb1;
    iov0 = iov1;
    goto L120;
L310:
    x1 = mgonxtv1_1.xx[1999];
    ya1 = mgoalin_(&x[ii], &x[ii + 1], &y[ii], &y[ii + 1], &x1);
    yb1 = mgonxtv1_1.yy[1999];
    if (ya1 > yb1) {
	goto L320;
    }
    mgooutp_(&x1, &yb1, ier);
    mgooutp_(&x1, &ya1, ier);
    px = x1;
    py = ya1;
    goto L330;
L380:
    ii = 1;
L320:
    px = x[ii];
    py = y[ii];
L330:
    if (ii == *n) {
	goto L400;
    }
    ++ii;
    mgooutp_(&x[ii], &y[ii], ier);
    mgoline_(&px, &py, &x[ii], &y[ii]);
    px = x[ii];
    py = y[ii];
    goto L330;
L300:
    if (jj == 2000) {
	goto L400;
    }
L340:
    x1 = x[*n];
    ya1 = y[*n];
    yb1 = mgoalin_(&mgonxtv1_1.xx[jj - 1], &mgonxtv1_1.xx[jj], &mgonxtv1_1.yy[
	    jj - 1], &mgonxtv1_1.yy[jj], &x1);
    if (ya1 <= yb1) {
	goto L350;
    }
    mgooutp_(&x1, &ya1, ier);
    mgooutp_(&x1, &yb1, ier);
    mgoline_(&xl, &yl, &x1, &ya1);
L350:
    if (jj == 2000) {
	goto L400;
    }
    ++jj;
    mgooutp_(&mgonxtv1_1.xx[jj - 1], &mgonxtv1_1.yy[jj - 1], ier);
    goto L350;
L360:
    jj = 0;
    goto L350;
L400:
    mgonxtv1_1.ll = 2000 - mgonxtv1_1.kk + 1;
    i__ = mgonxtv1_1.ll;
    i__1 = mgonxtv1_1.kk;
    for (j = 1; j <= i__1; ++j) {
	mgonxtv1_1.xx[i__ - 1] = mgonxtv1_1.xx[j - 1];
	mgonxtv1_1.yy[i__ - 1] = mgonxtv1_1.yy[j - 1];
	++i__;
/* L410: */
    }
    return 0;
L500:
    *ier = 1;
    return 0;
} /* mgonxtvu_ */

/* Subroutine */ int mgooutp_(x, y, ier)
real *x, *y;
integer *ier;
{
    /* System generated locals */
    real r__1, r__2;

    if (mgonxtv1_1.kk == 0) {
	goto L10;
    }
    if (mgonxtv1_1.kk == mgonxtv1_1.ll - 1) {
	goto L20;
    }
    if ((r__1 = mgonxtv1_1.xx[mgonxtv1_1.kk - 1] - *x, dabs(r__1)) + (r__2 = 
	    mgonxtv1_1.yy[mgonxtv1_1.kk - 1] - *y, dabs(r__2)) < (float).001) 
	    {
	return 0;
    }
L10:
    ++mgonxtv1_1.kk;
    mgonxtv1_1.xx[mgonxtv1_1.kk - 1] = *x;
    mgonxtv1_1.yy[mgonxtv1_1.kk - 1] = *y;
    return 0;
L20:
    *ier = 1;
    return 0;
} /* mgooutp_ */

doublereal mgoalin_(x0, x1, y0, y1, x)
real *x0, *x1, *y0, *y1, *x;
{
    /* System generated locals */
    real ret_val;

    if (*x0 == *x1) {
	goto L10;
    }
    ret_val = (*x - *x0) * (*y1 - *y0) / (*x1 - *x0) + *y0;
    return ret_val;
L10:
    if (*y1 > *y0) {
	goto L20;
    }
    ret_val = *y0;
    return ret_val;
L20:
    ret_val = *y1;
    return ret_val;
} /* mgoalin_ */

