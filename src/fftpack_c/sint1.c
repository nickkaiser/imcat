#include <math.h>
/* sint1.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int sint1_(n, war, was, xh, x, ifac)
integer *n;
real *war, *was, *xh, *x;
integer *ifac;
{
    /* Initialized data */

    static real sqrt3 = (float)1.73205080756888;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer modn, i__, k;
    static real xhold, t1, t2;
    extern /* Subroutine */ int rfftf1_();
    static integer kc, np1, ns2;

    /* Parameter adjustments */
    --ifac;
    --x;
    --xh;
    --was;
    --war;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xh[i__] = war[i__];
	war[i__] = x[i__];
/* L100: */
    }
    if ((i__1 = *n - 2) < 0) {
	goto L101;
    } else if (i__1 == 0) {
	goto L102;
    } else {
	goto L103;
    }
L101:
    xh[1] += xh[1];
    goto L106;
L102:
    xhold = sqrt3 * (xh[1] + xh[2]);
    xh[2] = sqrt3 * (xh[1] - xh[2]);
    xh[1] = xhold;
    goto L106;
L103:
    np1 = *n + 1;
    ns2 = *n / 2;
    x[1] = (float)0.;
    i__1 = ns2;
    for (k = 1; k <= i__1; ++k) {
	kc = np1 - k;
	t1 = xh[k] - xh[kc];
	t2 = was[k] * (xh[k] + xh[kc]);
	x[k + 1] = t1 + t2;
	x[kc + 1] = t2 - t1;
/* L104: */
    }
    modn = *n % 2;
    if (modn != 0) {
	x[ns2 + 2] = xh[ns2 + 1] * (float)4.;
    }
    rfftf1_(&np1, &x[1], &xh[1], &war[1], &ifac[1]);
    xh[1] = x[1] * (float).5;
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; i__ += 2) {
	xh[i__ - 1] = -x[i__];
	xh[i__] = xh[i__ - 2] + x[i__ - 1];
/* L105: */
    }
    if (modn != 0) {
	goto L106;
    }
    xh[*n] = -x[*n + 1];
L106:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = war[i__];
	war[i__] = xh[i__];
/* L107: */
    }
    return 0;
} /* sint1_ */

