#include <math.h>
/* radf3.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int radf3_(ido, l1, cc, ch, wa1, wa2)
integer *ido, *l1;
real *cc, *ch, *wa1, *wa2;
{
    /* Initialized data */

    static real taur = (float)-.5;
    static real taui = (float).866025403784439;

    /* System generated locals */
    integer ch_dim1, ch_offset, cc_dim1, cc_dim2, cc_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, ic;
    static real ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3;
    static integer idp2;

    /* Parameter adjustments */
    ch_dim1 = *ido;
    ch_offset = (ch_dim1 << 2) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_dim2 = *l1;
    cc_offset = cc_dim1 * (cc_dim2 + 1) + 1;
    cc -= cc_offset;
    --wa1;
    --wa2;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	cr2 = cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1] + cc[(k + cc_dim2 * 3) * 
		cc_dim1 + 1];
	ch[(k * 3 + 1) * ch_dim1 + 1] = cc[(k + cc_dim2) * cc_dim1 + 1] + cr2;
	ch[(k * 3 + 3) * ch_dim1 + 1] = taui * (cc[(k + cc_dim2 * 3) * 
		cc_dim1 + 1] - cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1]);
	ch[*ido + (k * 3 + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] + 
		taur * cr2;
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    dr2 = wa1[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1] 
		    + wa1[i__ - 1] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1];
	    di2 = wa1[i__ - 2] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1] - 
		    wa1[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * 
		    cc_dim1];
	    dr3 = wa2[i__ - 2] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1] + 
		    wa2[i__ - 1] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1];
	    di3 = wa2[i__ - 2] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1] - wa2[
		    i__ - 1] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1];
	    cr2 = dr2 + dr3;
	    ci2 = di2 + di3;
	    ch[i__ - 1 + (k * 3 + 1) * ch_dim1] = cc[i__ - 1 + (k + cc_dim2) *
		     cc_dim1] + cr2;
	    ch[i__ + (k * 3 + 1) * ch_dim1] = cc[i__ + (k + cc_dim2) * 
		    cc_dim1] + ci2;
	    tr2 = cc[i__ - 1 + (k + cc_dim2) * cc_dim1] + taur * cr2;
	    ti2 = cc[i__ + (k + cc_dim2) * cc_dim1] + taur * ci2;
	    tr3 = taui * (di2 - di3);
	    ti3 = taui * (dr3 - dr2);
	    ch[i__ - 1 + (k * 3 + 3) * ch_dim1] = tr2 + tr3;
	    ch[ic - 1 + (k * 3 + 2) * ch_dim1] = tr2 - tr3;
	    ch[i__ + (k * 3 + 3) * ch_dim1] = ti2 + ti3;
	    ch[ic + (k * 3 + 2) * ch_dim1] = ti3 - ti2;
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* radf3_ */

