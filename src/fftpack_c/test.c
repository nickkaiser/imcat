#include <math.h>
/* test.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c_n1 = -1;
static integer c__1 = 1;

/* Main program */ MAIN__()
{
    /* Initialized data */

    static integer nd[10] = { 120,54,49,32,4,3,2 };

    /* Format strings */
    static char fmt_1001[] = "(\0020N\002,i5,\002 RFFTF  \002,e10.3,\002 RFF\
TB  \002,e10.3,\002 RFFTFB \002,e10.3,\002 SINT   \002,e10.3,\002 SINTFB \
\002,e10.3,\002 COST   \002,e10.3/7x,\002 COSTFB \002,e10.3,\002 SINQF  \002\
,e10.3,\002 SINQB  \002,e10.3,\002 SINQFB \002,e10.3,\002 COSQF  \002,e10.3\
,\002 COSQB  \002,e10.3/7x,\002 COSQFB \002,e10.3,\002 DEZF   \002,e10.3,\
\002 DEZB   \002,e10.3,\002 DEZFB  \002,e10.3,\002 CFFTF  \002,e10.3,\002 CF\
FTB  \002,e10.3/7x,\002 CFFTFB \002,e10.3)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    double sqrt(), sin(), cos();
    integer pow_ii();
    double atan(), c_abs();
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer modn;
    static real rftb, rftf;
    extern /* Subroutine */ int cost_(), sint_();
    static real dezb1, dezf1, a[100], b[100];
    static integer i__, j;
    static real sqrt2;
    static integer k, n;
    extern /* Subroutine */ int cfftb_();
    static real w[2000], x[200], y[200], dezfb;
    extern /* Subroutine */ int cffti_(), cfftf_(), rfftb_();
    static real rftfb;
    extern /* Subroutine */ int rfftf_(), cosqb_(), rffti_(), cosqf_(), 
	    sinqb_(), cosqi_(), sinqf_(), costi_();
    static real azero;
    extern /* Subroutine */ int sinqi_(), sinti_();
    static real costt, sintt, ah[100], bh[100], cf, fn, dt, pi, dcfftb;
    static complex cx[200], cy[200];
    static real dcfftf, xh[200];
    static integer nz;
    static real cosqfb, costfb;
    extern /* Subroutine */ int ezfftb_();
    static real sinqfb;
    extern /* Subroutine */ int ezfftf_();
    static real sintfb;
    extern /* Subroutine */ int ezffti_();
    static real azeroh, cosqbt, cosqft, sinqbt, sinqft;
    static integer nm1, np1, ns2;
    static real arg, tfn, tpi;
    static integer nns;
    static real sum, arg1, arg2;
    static integer ns2m;
    static real sum1, sum2, dcfb;

    /* Fortran I/O blocks */
    static cilist io___57 = { 0, 6, 0, fmt_1001, 0 };



/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
*/

/*                       VERSION 4  APRIL 1985 */

/*                         A TEST DRIVER FOR */
/*          A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE FAST FOURIER */
/*           TRANSFORM OF PERIODIC AND OTHER SYMMETRIC SEQUENCES */

/*                              BY */

/*                       PAUL N SWARZTRAUBER */

/*       NATIONAL CENTER FOR ATMOSPHERIC RESEARCH  BOULDER,COLORADO 80307 
*/

/*        WHICH IS SPONSORED BY THE NATIONAL SCIENCE FOUNDATION */

/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
*/


/*             THIS PROGRAM TESTS THE PACKAGE OF FAST FOURIER */
/*     TRANSFORMS FOR BOTH COMPLEX AND REAL PERIODIC SEQUENCES AND */
/*     CERTIAN OTHER SYMMETRIC SEQUENCES THAT ARE LISTED BELOW. */

/*     1.   RFFTI     INITIALIZE  RFFTF AND RFFTB */
/*     2.   RFFTF     FORWARD TRANSFORM OF A REAL PERIODIC SEQUENCE */
/*     3.   RFFTB     BACKWARD TRANSFORM OF A REAL COEFFICIENT ARRAY */

/*     4.   EZFFTI    INITIALIZE EZFFTF AND EZFFTB */
/*     5.   EZFFTF    A SIMPLIFIED REAL PERIODIC FORWARD TRANSFORM */
/*     6.   EZFFTB    A SIMPLIFIED REAL PERIODIC BACKWARD TRANSFORM */

/*     7.   SINTI     INITIALIZE SINT */
/*     8.   SINT      SINE TRANSFORM OF A REAL ODD SEQUENCE */

/*     9.   COSTI     INITIALIZE COST */
/*     10.  COST      COSINE TRANSFORM OF A REAL EVEN SEQUENCE */

/*     11.  SINQI     INITIALIZE SINQF AND SINQB */
/*     12.  SINQF     FORWARD SINE TRANSFORM WITH ODD WAVE NUMBERS */
/*     13.  SINQB     UNNORMALIZED INVERSE OF SINQF */

/*     14.  COSQI     INITIALIZE COSQF AND COSQB */
/*     15.  COSQF     FORWARD COSINE TRANSFORM WITH ODD WAVE NUMBERS */
/*     16.  COSQB     UNNORMALIZED INVERSE OF COSQF */

/*     17.  CFFTI     INITIALIZE CFFTF AND CFFTB */
/*     18.  CFFTF     FORWARD TRANSFORM OF A COMPLEX PERIODIC SEQUENCE */
/*     19.  CFFTB     UNNORMALIZED INVERSE OF CFFTF */


    sqrt2 = sqrt((float)2.);
    nns = 7;
    i__1 = nns;
    for (nz = 1; nz <= i__1; ++nz) {
	n = nd[nz - 1];
	modn = n % 2;
	fn = (real) n;
	tfn = fn + fn;
	np1 = n + 1;
	nm1 = n - 1;
	i__2 = np1;
	for (j = 1; j <= i__2; ++j) {
	    x[j - 1] = sin((real) j * sqrt2);
	    y[j - 1] = x[j - 1];
	    xh[j - 1] = x[j - 1];
/* L101: */
	}

/*     TEST SUBROUTINES RFFTI,RFFTF AND RFFTB */

	rffti_(&n, w);
	pi = (float)3.14159265358979;
	dt = (pi + pi) / fn;
	ns2 = (n + 1) / 2;
	if (ns2 < 2) {
	    goto L104;
	}
	i__2 = ns2;
	for (k = 2; k <= i__2; ++k) {
	    sum1 = (float)0.;
	    sum2 = (float)0.;
	    arg = (real) (k - 1) * dt;
	    i__3 = n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		arg1 = (real) (i__ - 1) * arg;
		sum1 += x[i__ - 1] * cos(arg1);
		sum2 += x[i__ - 1] * sin(arg1);
/* L102: */
	    }
	    y[(k << 1) - 3] = sum1;
	    y[(k << 1) - 2] = -sum2;
/* L103: */
	}
L104:
	sum1 = (float)0.;
	sum2 = (float)0.;
	i__2 = nm1;
	for (i__ = 1; i__ <= i__2; i__ += 2) {
	    sum1 += x[i__ - 1];
	    sum2 += x[i__];
/* L105: */
	}
	if (modn == 1) {
	    sum1 += x[n - 1];
	}
	y[0] = sum1 + sum2;
	if (modn == 0) {
	    y[n - 1] = sum1 - sum2;
	}
	rfftf_(&n, x, w);
	rftf = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = rftf, r__3 = (r__1 = x[i__ - 1] - y[i__ - 1], dabs(r__1));
	    rftf = dmax(r__2,r__3);
	    x[i__ - 1] = xh[i__ - 1];
/* L106: */
	}
	rftf /= fn;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum = x[0] * (float).5;
	    arg = (real) (i__ - 1) * dt;
	    if (ns2 < 2) {
		goto L108;
	    }
	    i__3 = ns2;
	    for (k = 2; k <= i__3; ++k) {
		arg1 = (real) (k - 1) * arg;
		sum = sum + x[(k << 1) - 3] * cos(arg1) - x[(k << 1) - 2] * 
			sin(arg1);
/* L107: */
	    }
L108:
	    if (modn == 0) {
		i__3 = i__ - 1;
		sum += (real) pow_ii(&c_n1, &i__3) * (float).5 * x[n - 1];
	    }
	    y[i__ - 1] = sum + sum;
/* L109: */
	}
	rfftb_(&n, x, w);
	rftb = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = rftb, r__3 = (r__1 = x[i__ - 1] - y[i__ - 1], dabs(r__1));
	    rftb = dmax(r__2,r__3);
	    x[i__ - 1] = xh[i__ - 1];
	    y[i__ - 1] = xh[i__ - 1];
/* L110: */
	}
	rfftb_(&n, y, w);
	rfftf_(&n, y, w);
	cf = (float)1. / fn;
	rftfb = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = rftfb, r__3 = (r__1 = cf * y[i__ - 1] - x[i__ - 1], dabs(
		    r__1));
	    rftfb = dmax(r__2,r__3);
/* L111: */
	}

/*     TEST SUBROUTINES SINTI AND SINT */

	dt = pi / fn;
	i__2 = nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x[i__ - 1] = xh[i__ - 1];
/* L112: */
	}
	i__2 = nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ - 1] = (float)0.;
	    arg1 = (real) i__ * dt;
	    i__3 = nm1;
	    for (k = 1; k <= i__3; ++k) {
		y[i__ - 1] += x[k - 1] * sin((real) k * arg1);
/* L113: */
	    }
	    y[i__ - 1] += y[i__ - 1];
/* L114: */
	}
	sinti_(&nm1, w);
	sint_(&nm1, x, w);
	cf = (float).5 / fn;
	sintt = (float)0.;
	i__2 = nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = sintt, r__3 = (r__1 = x[i__ - 1] - y[i__ - 1], dabs(r__1));
	    sintt = dmax(r__2,r__3);
	    x[i__ - 1] = xh[i__ - 1];
	    y[i__ - 1] = x[i__ - 1];
/* L115: */
	}
	sintt = cf * sintt;
	sint_(&nm1, x, w);
	sint_(&nm1, x, w);
	sintfb = (float)0.;
	i__2 = nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = sintfb, r__3 = (r__1 = cf * x[i__ - 1] - y[i__ - 1], dabs(
		    r__1));
	    sintfb = dmax(r__2,r__3);
/* L116: */
	}

/*     TEST SUBROUTINES COSTI AND COST */

	i__2 = np1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x[i__ - 1] = xh[i__ - 1];
/* L117: */
	}
	i__2 = np1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + 1;
	    y[i__ - 1] = (x[0] + (real) pow_ii(&c_n1, &i__3) * x[n]) * (float)
		    .5;
	    arg = (real) (i__ - 1) * dt;
	    i__3 = n;
	    for (k = 2; k <= i__3; ++k) {
		y[i__ - 1] += x[k - 1] * cos((real) (k - 1) * arg);
/* L118: */
	    }
	    y[i__ - 1] += y[i__ - 1];
/* L119: */
	}
	costi_(&np1, w);
	cost_(&np1, x, w);
	costt = (float)0.;
	i__2 = np1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = costt, r__3 = (r__1 = x[i__ - 1] - y[i__ - 1], dabs(r__1));
	    costt = dmax(r__2,r__3);
	    x[i__ - 1] = xh[i__ - 1];
	    y[i__ - 1] = xh[i__ - 1];
/* L120: */
	}
	costt = cf * costt;
	cost_(&np1, x, w);
	cost_(&np1, x, w);
	costfb = (float)0.;
	i__2 = np1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = costfb, r__3 = (r__1 = cf * x[i__ - 1] - y[i__ - 1], dabs(
		    r__1));
	    costfb = dmax(r__2,r__3);
/* L121: */
	}

/*     TEST SUBROUTINES SINQI,SINQF AND SINQB */

	cf = (float).25 / fn;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ - 1] = xh[i__ - 1];
/* L122: */
	}
	dt = pi / (fn + fn);
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x[i__ - 1] = (float)0.;
	    arg = dt * (real) i__;
	    i__3 = n;
	    for (k = 1; k <= i__3; ++k) {
		x[i__ - 1] += y[k - 1] * sin((real) (k + k - 1) * arg);
/* L123: */
	    }
	    x[i__ - 1] *= (float)4.;
/* L124: */
	}
	sinqi_(&n, w);
	sinqb_(&n, y, w);
	sinqbt = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = sinqbt, r__3 = (r__1 = y[i__ - 1] - x[i__ - 1], dabs(r__1))
		    ;
	    sinqbt = dmax(r__2,r__3);
	    x[i__ - 1] = xh[i__ - 1];
/* L125: */
	}
	sinqbt = cf * sinqbt;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    arg = (real) (i__ + i__ - 1) * dt;
	    i__3 = i__ + 1;
	    y[i__ - 1] = (real) pow_ii(&c_n1, &i__3) * (float).5 * x[n - 1];
	    i__3 = nm1;
	    for (k = 1; k <= i__3; ++k) {
		y[i__ - 1] += x[k - 1] * sin((real) k * arg);
/* L126: */
	    }
	    y[i__ - 1] += y[i__ - 1];
/* L127: */
	}
	sinqf_(&n, x, w);
	sinqft = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = sinqft, r__3 = (r__1 = x[i__ - 1] - y[i__ - 1], dabs(r__1))
		    ;
	    sinqft = dmax(r__2,r__3);
	    y[i__ - 1] = xh[i__ - 1];
	    x[i__ - 1] = xh[i__ - 1];
/* L128: */
	}
	sinqf_(&n, y, w);
	sinqb_(&n, y, w);
	sinqfb = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = sinqfb, r__3 = (r__1 = cf * y[i__ - 1] - x[i__ - 1], dabs(
		    r__1));
	    sinqfb = dmax(r__2,r__3);
/* L129: */
	}

/*     TEST SUBROUTINES COSQI,COSQF AND COSQB */

	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ - 1] = xh[i__ - 1];
/* L130: */
	}
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x[i__ - 1] = (float)0.;
	    arg = (real) (i__ - 1) * dt;
	    i__3 = n;
	    for (k = 1; k <= i__3; ++k) {
		x[i__ - 1] += y[k - 1] * cos((real) (k + k - 1) * arg);
/* L131: */
	    }
	    x[i__ - 1] *= (float)4.;
/* L132: */
	}
	cosqi_(&n, w);
	cosqb_(&n, y, w);
	cosqbt = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = cosqbt, r__3 = (r__1 = x[i__ - 1] - y[i__ - 1], dabs(r__1))
		    ;
	    cosqbt = dmax(r__2,r__3);
	    x[i__ - 1] = xh[i__ - 1];
/* L133: */
	}
	cosqbt = cf * cosqbt;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ - 1] = x[0] * (float).5;
	    arg = (real) (i__ + i__ - 1) * dt;
	    i__3 = n;
	    for (k = 2; k <= i__3; ++k) {
		y[i__ - 1] += x[k - 1] * cos((real) (k - 1) * arg);
/* L134: */
	    }
	    y[i__ - 1] += y[i__ - 1];
/* L135: */
	}
	cosqf_(&n, x, w);
	cosqft = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = cosqft, r__3 = (r__1 = y[i__ - 1] - x[i__ - 1], dabs(r__1))
		    ;
	    cosqft = dmax(r__2,r__3);
	    x[i__ - 1] = xh[i__ - 1];
	    y[i__ - 1] = xh[i__ - 1];
/* L136: */
	}
	cosqft = cf * cosqft;
	cosqb_(&n, x, w);
	cosqf_(&n, x, w);
	cosqfb = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = cosqfb, r__3 = (r__1 = cf * x[i__ - 1] - y[i__ - 1], dabs(
		    r__1));
	    cosqfb = dmax(r__2,r__3);
/* L137: */
	}

/*     TEST PROGRAMS EZFFTI,EZFFTF,EZFFTB */

	ezffti_(&n, w);
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x[i__ - 1] = xh[i__ - 1];
/* L138: */
	}
	tpi = atan((float)1.) * (float)8.;
	dt = tpi / (real) n;
	ns2 = (n + 1) / 2;
	cf = (float)2. / (real) n;
	ns2m = ns2 - 1;
	if (ns2m <= 0) {
	    goto L141;
	}
	i__2 = ns2m;
	for (k = 1; k <= i__2; ++k) {
	    sum1 = (float)0.;
	    sum2 = (float)0.;
	    arg = (real) k * dt;
	    i__3 = n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		arg1 = (real) (i__ - 1) * arg;
		sum1 += x[i__ - 1] * cos(arg1);
		sum2 += x[i__ - 1] * sin(arg1);
/* L139: */
	    }
	    a[k - 1] = cf * sum1;
	    b[k - 1] = cf * sum2;
/* L140: */
	}
L141:
	nm1 = n - 1;
	sum1 = (float)0.;
	sum2 = (float)0.;
	i__2 = nm1;
	for (i__ = 1; i__ <= i__2; i__ += 2) {
	    sum1 += x[i__ - 1];
	    sum2 += x[i__];
/* L142: */
	}
	if (modn == 1) {
	    sum1 += x[n - 1];
	}
	azero = cf * (float).5 * (sum1 + sum2);
	if (modn == 0) {
	    a[ns2 - 1] = cf * (float).5 * (sum1 - sum2);
	}
	ezfftf_(&n, x, &azeroh, ah, bh, w);
	dezf1 = (r__1 = azeroh - azero, dabs(r__1));
	if (modn == 0) {
/* Computing MAX */
	    r__2 = dezf1, r__3 = (r__1 = a[ns2 - 1] - ah[ns2 - 1], dabs(r__1))
		    ;
	    dezf1 = dmax(r__2,r__3);
	}
	if (ns2m <= 0) {
	    goto L144;
	}
	i__2 = ns2m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__3 = dezf1, r__4 = (r__1 = ah[i__ - 1] - a[i__ - 1], dabs(r__1))
		    , r__3 = max(r__3,r__4), r__4 = (r__2 = bh[i__ - 1] - b[
		    i__ - 1], dabs(r__2));
	    dezf1 = dmax(r__3,r__4);
/* L143: */
	}
L144:
	ns2 = n / 2;
	if (modn == 0) {
	    b[ns2 - 1] = (float)0.;
	}
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum = azero;
	    arg1 = (real) (i__ - 1) * dt;
	    i__3 = ns2;
	    for (k = 1; k <= i__3; ++k) {
		arg2 = (real) k * arg1;
		sum = sum + a[k - 1] * cos(arg2) + b[k - 1] * sin(arg2);
/* L145: */
	    }
	    x[i__ - 1] = sum;
/* L146: */
	}
	ezfftb_(&n, y, &azero, a, b, w);
	dezb1 = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = dezb1, r__3 = (r__1 = x[i__ - 1] - y[i__ - 1], dabs(r__1));
	    dezb1 = dmax(r__2,r__3);
	    x[i__ - 1] = xh[i__ - 1];
/* L147: */
	}
	ezfftf_(&n, x, &azero, a, b, w);
	ezfftb_(&n, y, &azero, a, b, w);
	dezfb = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = dezfb, r__3 = (r__1 = x[i__ - 1] - y[i__ - 1], dabs(r__1));
	    dezfb = dmax(r__2,r__3);
/* L148: */
	}

/*     TEST  CFFTI,CFFTF,CFFTB */

	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ - 1;
	    r__1 = cos(sqrt2 * (real) i__);
	    r__2 = sin(sqrt2 * (real) (i__ * i__));
	    q__1.r = r__1, q__1.i = r__2;
	    cx[i__3].r = q__1.r, cx[i__3].i = q__1.i;
/* L149: */
	}
	dt = (pi + pi) / fn;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    arg1 = -((real) (i__ - 1)) * dt;
	    i__3 = i__ - 1;
	    cy[i__3].r = (float)0., cy[i__3].i = (float)0.;
	    i__3 = n;
	    for (k = 1; k <= i__3; ++k) {
		arg2 = (real) (k - 1) * arg1;
		i__4 = i__ - 1;
		i__5 = i__ - 1;
		r__1 = cos(arg2);
		r__2 = sin(arg2);
		q__3.r = r__1, q__3.i = r__2;
		i__6 = k - 1;
		q__2.r = q__3.r * cx[i__6].r - q__3.i * cx[i__6].i, q__2.i = 
			q__3.r * cx[i__6].i + q__3.i * cx[i__6].r;
		q__1.r = cy[i__5].r + q__2.r, q__1.i = cy[i__5].i + q__2.i;
		cy[i__4].r = q__1.r, cy[i__4].i = q__1.i;
/* L150: */
	    }
/* L151: */
	}
	cffti_(&n, w);
	cfftf_(&n, cx, w);
	dcfftf = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    i__3 = i__ - 1;
	    i__4 = i__ - 1;
	    q__1.r = cx[i__3].r - cy[i__4].r, q__1.i = cx[i__3].i - cy[i__4]
		    .i;
	    r__1 = dcfftf, r__2 = c_abs(&q__1);
	    dcfftf = dmax(r__1,r__2);
	    i__3 = i__ - 1;
	    i__4 = i__ - 1;
	    q__1.r = cx[i__4].r / fn, q__1.i = cx[i__4].i / fn;
	    cx[i__3].r = q__1.r, cx[i__3].i = q__1.i;
/* L152: */
	}
	dcfftf /= fn;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    arg1 = (real) (i__ - 1) * dt;
	    i__3 = i__ - 1;
	    cy[i__3].r = (float)0., cy[i__3].i = (float)0.;
	    i__3 = n;
	    for (k = 1; k <= i__3; ++k) {
		arg2 = (real) (k - 1) * arg1;
		i__4 = i__ - 1;
		i__5 = i__ - 1;
		r__1 = cos(arg2);
		r__2 = sin(arg2);
		q__3.r = r__1, q__3.i = r__2;
		i__6 = k - 1;
		q__2.r = q__3.r * cx[i__6].r - q__3.i * cx[i__6].i, q__2.i = 
			q__3.r * cx[i__6].i + q__3.i * cx[i__6].r;
		q__1.r = cy[i__5].r + q__2.r, q__1.i = cy[i__5].i + q__2.i;
		cy[i__4].r = q__1.r, cy[i__4].i = q__1.i;
/* L153: */
	    }
/* L154: */
	}
	cfftb_(&n, cx, w);
	dcfftb = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    i__3 = i__ - 1;
	    i__4 = i__ - 1;
	    q__1.r = cx[i__3].r - cy[i__4].r, q__1.i = cx[i__3].i - cy[i__4]
		    .i;
	    r__1 = dcfftb, r__2 = c_abs(&q__1);
	    dcfftb = dmax(r__1,r__2);
	    i__3 = i__ - 1;
	    i__4 = i__ - 1;
	    cx[i__3].r = cy[i__4].r, cx[i__3].i = cy[i__4].i;
/* L155: */
	}
	cf = (float)1. / fn;
	cfftf_(&n, cx, w);
	cfftb_(&n, cx, w);
	dcfb = (float)0.;
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    i__3 = i__ - 1;
	    q__2.r = cf * cx[i__3].r, q__2.i = cf * cx[i__3].i;
	    i__4 = i__ - 1;
	    q__1.r = q__2.r - cy[i__4].r, q__1.i = q__2.i - cy[i__4].i;
	    r__1 = dcfb, r__2 = c_abs(&q__1);
	    dcfb = dmax(r__1,r__2);
/* L156: */
	}
	s_wsfe(&io___57);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&rftf, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rftb, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rftfb, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&sintt, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&sintfb, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&costt, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&costfb, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&sinqft, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&sinqbt, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&sinqfb, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cosqft, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cosqbt, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cosqfb, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&dezf1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&dezb1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&dezfb, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&dcfftf, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&dcfftb, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&dcfb, (ftnlen)sizeof(real));
	e_wsfe();
/* L157: */
    }




} /* MAIN__ */

/* Main program alias */ int tstfft_ () { MAIN__ (); }
