#include <math.h>
/* ezfftf.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int ezfftf_(n, r__, azero, a, b, wsave)
integer *n;
real *r__, *azero, *a, *b, *wsave;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int rfftf_();
    static real cf;
    static integer ns2;
    static real cfm;
    static integer ns2m;


/*                       VERSION 3  JUNE 1979 */

    /* Parameter adjustments */
    --wsave;
    --b;
    --a;
    --r__;

    /* Function Body */
    if ((i__1 = *n - 2) < 0) {
	goto L101;
    } else if (i__1 == 0) {
	goto L102;
    } else {
	goto L103;
    }
L101:
    *azero = r__[1];
    return 0;
L102:
    *azero = (r__[1] + r__[2]) * (float).5;
    a[1] = (r__[1] - r__[2]) * (float).5;
    return 0;
L103:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wsave[i__] = r__[i__];
/* L104: */
    }
    rfftf_(n, &wsave[1], &wsave[*n + 1]);
    cf = (float)2. / (real) (*n);
    cfm = -cf;
    *azero = cf * (float).5 * wsave[1];
    ns2 = (*n + 1) / 2;
    ns2m = ns2 - 1;
    i__1 = ns2m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__] = cf * wsave[i__ * 2];
	b[i__] = cfm * wsave[(i__ << 1) + 1];
/* L105: */
    }
    if (*n % 2 == 1) {
	return 0;
    }
    a[ns2] = cf * (float).5 * wsave[*n];
    b[ns2] = (float)0.;
    return 0;
} /* ezfftf_ */

