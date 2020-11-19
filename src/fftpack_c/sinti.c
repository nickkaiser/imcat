#include <math.h>
/* sinti.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int sinti_(n, wsave)
integer *n;
real *wsave;
{
    /* Initialized data */

    static real pi = (float)3.14159265358979;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin();

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int rffti_();
    static real dt;
    static integer np1, ns2;

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n <= 1) {
	return 0;
    }
    ns2 = *n / 2;
    np1 = *n + 1;
    dt = pi / (real) np1;
    i__1 = ns2;
    for (k = 1; k <= i__1; ++k) {
	wsave[k] = sin(k * dt) * (float)2.;
/* L101: */
    }
    rffti_(&np1, &wsave[ns2 + 1]);
    return 0;
} /* sinti_ */

