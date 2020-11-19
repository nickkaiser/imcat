#include <math.h>
/* sint.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int sint_(n, x, wsave)
integer *n;
real *x, *wsave;
{
    extern /* Subroutine */ int sint1_();
    static integer np1, iw1, iw2, iw3;

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    np1 = *n + 1;
    iw1 = *n / 2 + 1;
    iw2 = iw1 + np1;
    iw3 = iw2 + np1;
    sint1_(n, &x[1], &wsave[1], &wsave[iw1], &wsave[iw2], &wsave[iw3]);
    return 0;
} /* sint_ */

