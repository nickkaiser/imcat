#include <math.h>
/* ezfftb.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int ezfftb_(n, r__, azero, a, b, wsave)
integer *n;
real *r__, *azero, *a, *b, *wsave;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int rfftb_();
    static integer ns2;

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
    r__[1] = *azero;
    return 0;
L102:
    r__[1] = *azero + a[1];
    r__[2] = *azero - a[1];
    return 0;
L103:
    ns2 = (*n - 1) / 2;
    i__1 = ns2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__[i__ * 2] = a[i__] * (float).5;
	r__[(i__ << 1) + 1] = b[i__] * (float)-.5;
/* L104: */
    }
    r__[1] = *azero;
    if (*n % 2 == 0) {
	r__[*n] = a[ns2 + 1];
    }
    rfftb_(n, &r__[1], &wsave[*n + 1]);
    return 0;
} /* ezfftb_ */

