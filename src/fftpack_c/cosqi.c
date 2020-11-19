#include <math.h>
/* cosqi.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int cosqi_(n, wsave)
integer *n;
real *wsave;
{
    /* Initialized data */

    static real pih = (float)1.57079632679491;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos();

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int rffti_();
    static real fk, dt;

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    dt = pih / (real) (*n);
    fk = (float)0.;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	fk += (float)1.;
	wsave[k] = cos(fk * dt);
/* L101: */
    }
    rffti_(n, &wsave[*n + 1]);
    return 0;
} /* cosqi_ */

