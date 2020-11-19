#include <math.h>
/* cfftf1.f -- translated by f2c (version 19951025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int cfftf1_(n, c__, ch, wa, ifac)
integer *n;
real *c__, *ch, *wa;
integer *ifac;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer idot, i__;
    extern /* Subroutine */ int passf_();
    static integer k1, l1, l2, n2;
    extern /* Subroutine */ int passf2_(), passf3_(), passf4_(), passf5_();
    static integer na, nf, ip, iw, ix2, ix3, ix4, nac, ido, idl1;

    /* Parameter adjustments */
    --ifac;
    --wa;
    --ch;
    --c__;

    /* Function Body */
    nf = ifac[2];
    na = 0;
    l1 = 1;
    iw = 1;
    i__1 = nf;
    for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	l2 = ip * l1;
	ido = *n / l2;
	idot = ido + ido;
	idl1 = idot * l1;
	if (ip != 4) {
	    goto L103;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	if (na != 0) {
	    goto L101;
	}
	passf4_(&idot, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L102;
L101:
	passf4_(&idot, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3]);
L102:
	na = 1 - na;
	goto L115;
L103:
	if (ip != 2) {
	    goto L106;
	}
	if (na != 0) {
	    goto L104;
	}
	passf2_(&idot, &l1, &c__[1], &ch[1], &wa[iw]);
	goto L105;
L104:
	passf2_(&idot, &l1, &ch[1], &c__[1], &wa[iw]);
L105:
	na = 1 - na;
	goto L115;
L106:
	if (ip != 3) {
	    goto L109;
	}
	ix2 = iw + idot;
	if (na != 0) {
	    goto L107;
	}
	passf3_(&idot, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L108;
L107:
	passf3_(&idot, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2]);
L108:
	na = 1 - na;
	goto L115;
L109:
	if (ip != 5) {
	    goto L112;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	ix4 = ix3 + idot;
	if (na != 0) {
	    goto L110;
	}
	passf5_(&idot, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
	goto L111;
L110:
	passf5_(&idot, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
L111:
	na = 1 - na;
	goto L115;
L112:
	if (na != 0) {
	    goto L113;
	}
	passf_(&nac, &idot, &ip, &l1, &idl1, &c__[1], &c__[1], &c__[1], &ch[1]
		, &ch[1], &wa[iw]);
	goto L114;
L113:
	passf_(&nac, &idot, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c__[1], 
		&c__[1], &wa[iw]);
L114:
	if (nac != 0) {
	    na = 1 - na;
	}
L115:
	l1 = l2;
	iw += (ip - 1) * idot;
/* L116: */
    }
    if (na == 0) {
	return 0;
    }
    n2 = *n + *n;
    i__1 = n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = ch[i__];
/* L117: */
    }
    return 0;
} /* cfftf1_ */

