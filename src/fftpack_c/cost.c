#include <math.h>
/* cost.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int cost_(n, x, wsave)
integer *n;
real *x, *wsave;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer modn, i__, k;
    extern /* Subroutine */ int rfftf_();
    static real c1, t1, t2;
    static integer kc;
    static real xi;
    static integer nm1, np1;
    static real x1h;
    static integer ns2;
    static real tx2, x1p3, xim2;

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    nm1 = *n - 1;
    np1 = *n + 1;
    ns2 = *n / 2;
    if ((i__1 = *n - 2) < 0) {
	goto L106;
    } else if (i__1 == 0) {
	goto L101;
    } else {
	goto L102;
    }
L101:
    x1h = x[1] + x[2];
    x[2] = x[1] - x[2];
    x[1] = x1h;
    return 0;
L102:
    if (*n > 3) {
	goto L103;
    }
    x1p3 = x[1] + x[3];
    tx2 = x[2] + x[2];
    x[2] = x[1] - x[3];
    x[1] = x1p3 + tx2;
    x[3] = x1p3 - tx2;
    return 0;
L103:
    c1 = x[1] - x[*n];
    x[1] += x[*n];
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np1 - k;
	t1 = x[k] + x[kc];
	t2 = x[k] - x[kc];
	c1 += wsave[kc] * t2;
	t2 = wsave[k] * t2;
	x[k] = t1 - t2;
	x[kc] = t1 + t2;
/* L104: */
    }
    modn = *n % 2;
    if (modn != 0) {
	x[ns2 + 1] += x[ns2 + 1];
    }
    rfftf_(&nm1, &x[1], &wsave[*n + 1]);
    xim2 = x[2];
    x[2] = c1;
    i__1 = *n;
    for (i__ = 4; i__ <= i__1; i__ += 2) {
	xi = x[i__];
	x[i__] = x[i__ - 2] - x[i__ - 1];
	x[i__ - 1] = xim2;
	xim2 = xi;
/* L105: */
    }
    if (modn != 0) {
	x[*n] = xim2;
    }
L106:
    return 0;
} /* cost_ */

