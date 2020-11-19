#include <math.h>
/* cosqb.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int cosqb_(n, x, wsave)
integer *n;
real *x, *wsave;
{
    /* Initialized data */

    static real tsqrt2 = (float)2.82842712474619;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real x1;
    extern /* Subroutine */ int cosqb1_();

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    if ((i__1 = *n - 2) < 0) {
	goto L101;
    } else if (i__1 == 0) {
	goto L102;
    } else {
	goto L103;
    }
L101:
    x[1] *= (float)4.;
    return 0;
L102:
    x1 = (x[1] + x[2]) * (float)4.;
    x[2] = tsqrt2 * (x[1] - x[2]);
    x[1] = x1;
    return 0;
L103:
    cosqb1_(n, &x[1], &wsave[1], &wsave[*n + 1]);
    return 0;
} /* cosqb_ */

