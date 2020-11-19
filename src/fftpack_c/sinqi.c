#include <math.h>
/* sinqi.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int sinqi_(n, wsave)
integer *n;
real *wsave;
{
    extern /* Subroutine */ int cosqi_();

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    cosqi_(n, &wsave[1]);
    return 0;
} /* sinqi_ */

