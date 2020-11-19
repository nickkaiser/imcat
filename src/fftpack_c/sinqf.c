#include <math.h>
/* sinqf.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int sinqf_(n, x, wsave)
integer *n;
real *x, *wsave;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int cosqf_();
    static real xhold;
    static integer kc, ns2;

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    ns2 = *n / 2;
    i__1 = ns2;
    for (k = 1; k <= i__1; ++k) {
	kc = *n - k;
	xhold = x[k];
	x[k] = x[kc + 1];
	x[kc + 1] = xhold;
/* L101: */
    }
    cosqf_(n, &x[1], &wsave[1]);
    i__1 = *n;
    for (k = 2; k <= i__1; k += 2) {
	x[k] = -x[k];
/* L102: */
    }
    return 0;
} /* sinqf_ */

