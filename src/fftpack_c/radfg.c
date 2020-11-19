#include <math.h>
/* radfg.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int radfg_(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
integer *ido, *ip, *l1, *idl1;
real *cc, *c1, *c2, *ch, *ch2, *wa;
{
    /* Initialized data */

    static real tpi = (float)6.28318530717959;

    /* System generated locals */
    integer ch_dim1, ch_dim2, ch_offset, cc_dim1, cc_dim2, cc_offset, c1_dim1,
	     c1_dim2, c1_offset, c2_dim1, c2_offset, ch2_dim1, ch2_offset, 
	    i__1, i__2, i__3;

    /* Builtin functions */
    double cos(), sin();

    /* Local variables */
    static integer idij, ipph, i__, j, k, l, j2, ic, jc, lc, ik, is;
    static real dc2, ai1, ai2, ar1, ar2, ds2;
    static integer nbd;
    static real dcp, arg, dsp, ar1h, ar2h;
    static integer idp2, ipp2;

    /* Parameter adjustments */
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    c1_dim1 = *ido;
    c1_dim2 = *l1;
    c1_offset = c1_dim1 * (c1_dim2 + 1) + 1;
    c1 -= c1_offset;
    cc_dim1 = *ido;
    cc_dim2 = *ip;
    cc_offset = cc_dim1 * (cc_dim2 + 1) + 1;
    cc -= cc_offset;
    ch2_dim1 = *idl1;
    ch2_offset = ch2_dim1 + 1;
    ch2 -= ch2_offset;
    c2_dim1 = *idl1;
    c2_offset = c2_dim1 + 1;
    c2 -= c2_offset;
    --wa;

    /* Function Body */
    arg = tpi / (real) (*ip);
    dcp = cos(arg);
    dsp = sin(arg);
    ipph = (*ip + 1) / 2;
    ipp2 = *ip + 2;
    idp2 = *ido + 2;
    nbd = (*ido - 1) / 2;
    if (*ido == 1) {
	goto L119;
    }
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	ch2[ik + ch2_dim1] = c2[ik + c2_dim1];
/* L101: */
    }
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[(k + j * ch_dim2) * ch_dim1 + 1] = c1[(k + j * c1_dim2) * 
		    c1_dim1 + 1];
/* L102: */
	}
/* L103: */
    }
    if (nbd > *l1) {
	goto L107;
    }
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	idij = is;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[
			i__ - 1 + (k + j * c1_dim2) * c1_dim1] + wa[idij] * 
			c1[i__ + (k + j * c1_dim2) * c1_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[i__ 
			+ (k + j * c1_dim2) * c1_dim1] - wa[idij] * c1[i__ - 
			1 + (k + j * c1_dim2) * c1_dim1];
/* L104: */
	    }
/* L105: */
	}
/* L106: */
    }
    goto L111;
L107:
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = is;
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		idij += 2;
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[
			i__ - 1 + (k + j * c1_dim2) * c1_dim1] + wa[idij] * 
			c1[i__ + (k + j * c1_dim2) * c1_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[i__ 
			+ (k + j * c1_dim2) * c1_dim1] - wa[idij] * c1[i__ - 
			1 + (k + j * c1_dim2) * c1_dim1];
/* L108: */
	    }
/* L109: */
	}
/* L110: */
    }
L111:
    if (nbd < *l1) {
	goto L115;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = ch[i__ - 1 + (k + 
			j * ch_dim2) * ch_dim1] + ch[i__ - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		c1[i__ - 1 + (k + jc * c1_dim2) * c1_dim1] = ch[i__ + (k + j *
			 ch_dim2) * ch_dim1] - ch[i__ + (k + jc * ch_dim2) * 
			ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = ch[i__ + (k + j * 
			ch_dim2) * ch_dim1] + ch[i__ + (k + jc * ch_dim2) * 
			ch_dim1];
		c1[i__ + (k + jc * c1_dim2) * c1_dim1] = ch[i__ - 1 + (k + jc 
			* ch_dim2) * ch_dim1] - ch[i__ - 1 + (k + j * ch_dim2)
			 * ch_dim1];
/* L112: */
	    }
/* L113: */
	}
/* L114: */
    }
    goto L121;
L115:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = ch[i__ - 1 + (k + 
			j * ch_dim2) * ch_dim1] + ch[i__ - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		c1[i__ - 1 + (k + jc * c1_dim2) * c1_dim1] = ch[i__ + (k + j *
			 ch_dim2) * ch_dim1] - ch[i__ + (k + jc * ch_dim2) * 
			ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = ch[i__ + (k + j * 
			ch_dim2) * ch_dim1] + ch[i__ + (k + jc * ch_dim2) * 
			ch_dim1];
		c1[i__ + (k + jc * c1_dim2) * c1_dim1] = ch[i__ - 1 + (k + jc 
			* ch_dim2) * ch_dim1] - ch[i__ - 1 + (k + j * ch_dim2)
			 * ch_dim1];
/* L116: */
	    }
/* L117: */
	}
/* L118: */
    }
    goto L121;
L119:
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	c2[ik + c2_dim1] = ch2[ik + ch2_dim1];
/* L120: */
    }
L121:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    c1[(k + j * c1_dim2) * c1_dim1 + 1] = ch[(k + j * ch_dim2) * 
		    ch_dim1 + 1] + ch[(k + jc * ch_dim2) * ch_dim1 + 1];
	    c1[(k + jc * c1_dim2) * c1_dim1 + 1] = ch[(k + jc * ch_dim2) * 
		    ch_dim1 + 1] - ch[(k + j * ch_dim2) * ch_dim1 + 1];
/* L122: */
	}
/* L123: */
    }

    ar1 = (float)1.;
    ai1 = (float)0.;
    i__1 = ipph;
    for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	ar1h = dcp * ar1 - dsp * ai1;
	ai1 = dcp * ai1 + dsp * ar1;
	ar1 = ar1h;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + l * ch2_dim1] = c2[ik + c2_dim1] + ar1 * c2[ik + (
		    c2_dim1 << 1)];
	    ch2[ik + lc * ch2_dim1] = ai1 * c2[ik + *ip * c2_dim1];
/* L124: */
	}
	dc2 = ar1;
	ds2 = ai1;
	ar2 = ar1;
	ai2 = ai1;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    ar2h = dc2 * ar2 - ds2 * ai2;
	    ai2 = dc2 * ai2 + ds2 * ar2;
	    ar2 = ar2h;
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		ch2[ik + l * ch2_dim1] += ar2 * c2[ik + j * c2_dim1];
		ch2[ik + lc * ch2_dim1] += ai2 * c2[ik + jc * c2_dim1];
/* L125: */
	    }
/* L126: */
	}
/* L127: */
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + ch2_dim1] += c2[ik + j * c2_dim1];
/* L128: */
	}
/* L129: */
    }

    if (*ido < *l1) {
	goto L132;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    cc[i__ + (k * cc_dim2 + 1) * cc_dim1] = ch[i__ + (k + ch_dim2) * 
		    ch_dim1];
/* L130: */
	}
/* L131: */
    }
    goto L135;
L132:
    i__1 = *ido;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    cc[i__ + (k * cc_dim2 + 1) * cc_dim1] = ch[i__ + (k + ch_dim2) * 
		    ch_dim1];
/* L133: */
	}
/* L134: */
    }
L135:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    cc[*ido + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[(k + j * ch_dim2)
		     * ch_dim1 + 1];
	    cc[(j2 - 1 + k * cc_dim2) * cc_dim1 + 1] = ch[(k + jc * ch_dim2) *
		     ch_dim1 + 1];
/* L136: */
	}
/* L137: */
    }
    if (*ido == 1) {
	return 0;
    }
    if (nbd < *l1) {
	goto L141;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		ic = idp2 - i__;
		cc[i__ - 1 + (j2 - 1 + k * cc_dim2) * cc_dim1] = ch[i__ - 1 + 
			(k + j * ch_dim2) * ch_dim1] + ch[i__ - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		cc[ic - 1 + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[i__ - 1 + (
			k + j * ch_dim2) * ch_dim1] - ch[i__ - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		cc[i__ + (j2 - 1 + k * cc_dim2) * cc_dim1] = ch[i__ + (k + j *
			 ch_dim2) * ch_dim1] + ch[i__ + (k + jc * ch_dim2) * 
			ch_dim1];
		cc[ic + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[i__ + (k + jc *
			 ch_dim2) * ch_dim1] - ch[i__ + (k + j * ch_dim2) * 
			ch_dim1];
/* L138: */
	    }
/* L139: */
	}
/* L140: */
    }
    return 0;
L141:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		cc[i__ - 1 + (j2 - 1 + k * cc_dim2) * cc_dim1] = ch[i__ - 1 + 
			(k + j * ch_dim2) * ch_dim1] + ch[i__ - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		cc[ic - 1 + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[i__ - 1 + (
			k + j * ch_dim2) * ch_dim1] - ch[i__ - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		cc[i__ + (j2 - 1 + k * cc_dim2) * cc_dim1] = ch[i__ + (k + j *
			 ch_dim2) * ch_dim1] + ch[i__ + (k + jc * ch_dim2) * 
			ch_dim1];
		cc[ic + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[i__ + (k + jc *
			 ch_dim2) * ch_dim1] - ch[i__ + (k + j * ch_dim2) * 
			ch_dim1];
/* L142: */
	    }
/* L143: */
	}
/* L144: */
    }
    return 0;
} /* radfg_ */

