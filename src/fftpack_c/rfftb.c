#include <math.h>
/* rfftb.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int rfftb_(n, r__, wsave)
integer *n;
real *r__, *wsave;
{
    extern /* Subroutine */ int rfftb1_();

    /* Parameter adjustments */
    --wsave;
    --r__;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    rfftb1_(n, &r__[1], &wsave[1], &wsave[*n + 1], &wsave[(*n << 1) + 1]);
    return 0;
} /* rfftb_ */

