#include <math.h>
/* cfftb.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int cfftb_(n, c__, wsave)
integer *n;
real *c__, *wsave;
{
    extern /* Subroutine */ int cfftb1_();
    static integer iw1, iw2;

    /* Parameter adjustments */
    --wsave;
    --c__;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    iw1 = *n + *n + 1;
    iw2 = iw1 + *n + *n;
    cfftb1_(n, &c__[1], &wsave[1], &wsave[iw1], &wsave[iw2]);
    return 0;
} /* cfftb_ */

