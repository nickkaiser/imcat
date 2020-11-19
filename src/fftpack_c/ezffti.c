#include <math.h>
/* ezffti.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int ezffti_(n, wsave)
integer *n;
real *wsave;
{
    extern /* Subroutine */ int ezfft1_();

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    ezfft1_(n, &wsave[(*n << 1) + 1], &wsave[*n * 3 + 1]);
    return 0;
} /* ezffti_ */

