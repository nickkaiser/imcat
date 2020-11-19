#include <math.h>
/* cffti.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int cffti_(n, wsave)
integer *n;
real *wsave;
{
    extern /* Subroutine */ int cffti1_();
    static integer iw1, iw2;

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    iw1 = *n + *n + 1;
    iw2 = iw1 + *n + *n;
    cffti1_(n, &wsave[iw1], &wsave[iw2]);
    return 0;
} /* cffti_ */

