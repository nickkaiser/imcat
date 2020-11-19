#include <math.h>
/* cosqf.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int cosqf_(n, x, wsave)
integer *n;
real *x, *wsave;
{
    /* Initialized data */

    static real sqrt2 = (float)1.4142135623731;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real tsqx;
    extern /* Subroutine */ int cosqf1_();

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    if ((i__1 = *n - 2) < 0) {
	goto L102;
    } else if (i__1 == 0) {
	goto L101;
    } else {
	goto L103;
    }
L101:
    tsqx = sqrt2 * x[2];
    x[2] = x[1] - tsqx;
    x[1] += tsqx;
L102:
    return 0;
L103:
    cosqf1_(n, &x[1], &wsave[1], &wsave[*n + 1]);
    return 0;
} /* cosqf_ */

