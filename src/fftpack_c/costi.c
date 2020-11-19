#include <math.h>
/* costi.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int costi_(n, wsave)
integer *n;
real *wsave;
{
    /* Initialized data */

    static real pi = (float)3.14159265358979;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(), cos();

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int rffti_();
    static integer kc;
    static real fk, dt;
    static integer nm1, np1, ns2;

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n <= 3) {
	return 0;
    }
    nm1 = *n - 1;
    np1 = *n + 1;
    ns2 = *n / 2;
    dt = pi / (real) nm1;
    fk = (float)0.;
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np1 - k;
	fk += (float)1.;
	wsave[k] = sin(fk * dt) * (float)2.;
	wsave[kc] = cos(fk * dt) * (float)2.;
/* L101: */
    }
    rffti_(&nm1, &wsave[*n + 1]);
    return 0;
} /* costi_ */

