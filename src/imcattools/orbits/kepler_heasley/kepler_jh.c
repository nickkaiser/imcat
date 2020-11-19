/* kepler_jh.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct consts_1_ {
    doublereal pi, twopi, xmu, dummy[4];
};
struct consts_2_ {
    doublereal pi, twopi, xmu, xk, eps, c__, rearth;
};

#define consts_1 (*(struct consts_1_ *) &consts_)
#define consts_2 (*(struct consts_2_ *) &consts_)

/* Initialized data */

struct {
    doublereal e_1[7];
    } consts_ = { 3.1415926535897931, 6.2831853071795862, 2.959122370203e-4, 
	    .01720209895, .40927971, 173.1428, 4.25493e-5 };


/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;


/* 	--------------------------------------------------- */
/* 	NK: changed to a integer function */
/* 	subroutine elmnts(elem,emu,t,r,rd) */
integer ilmnts_(doublereal *elem, doublereal *emu, doublereal *t, doublereal *
	r__, doublereal *rd)
{
    /* System generated locals */
    integer ret_val;
    doublereal d__1;

    /* Builtin functions */
    double acos(doublereal), sin(doublereal), atan2(doublereal, doublereal), 
	    sqrt(doublereal);

    /* Local variables */
    extern doublereal vlen_(doublereal *), vdot_(doublereal *, doublereal *);
    static doublereal c__, h__[3], s, u, v, ea, hm, cv, rm, sv, cea, sea, rdm,
	     hsi, per, csu;
    extern /* Subroutine */ int vxp_(doublereal *, doublereal *, doublereal *)
	    ;


/* 	given r and rdot vectors determine the orbital elements */


/* 	compute components of h */

/* 	h(1) = r(2)*rd(3) - r(3)*rd(2) */
/* 	h(2) = r(3)*rd(1) - r(1)*rd(3) */
/* 	h(3) = r(1)*rd(2) - r(2)*rd(1) */
/* 	hm = dsqrt(vdot(h,h)) */
    vxp_(r__, rd, h__);
    hm = vlen_(h__);
    elem[2] = acos(h__[2] / hm);
    hsi = hm * sin(elem[2]);
    c__ = -h__[1] / hsi;
    s = h__[0] / hsi;
    elem[3] = atan2(s, c__);
    if (elem[3] < 0.f) {
	elem[3] += consts_1.twopi;
    }
/* 	rm = dsqrt(vdot(r,r)) */
/* 	rdm = dsqrt(vdot(rd,rd)) */
    rm = vlen_(r__);
    rdm = vlen_(rd);
    elem[0] = *emu * rm / (*emu * 2.f - rm * rdm * rdm);
    if (elem[0] < 0.) {
/* 		write(*,2) */
/* 2		format(' orbit is hyperboic') */
/* 		stop */
	ret_val = 1;
	return ret_val;
    }
    elem[1] = sqrt(1.f - hm * hm / (*emu * elem[0]));
    csu = (c__ * r__[0] + s * r__[1]) / rm;
    u = acos(csu);
    if (r__[2] < 0.f) {
	u = consts_1.twopi - u;
    }
    cv = (hm * hm / (*emu * rm) - 1.f) / elem[1];
    sv = hm * vdot_(r__, rd) / (elem[1] * *emu * rm);
    v = atan2(sv, cv);
    if (v < 0.) {
	v += consts_1.twopi;
    }
    elem[4] = u - v;
    if (elem[4] < 0.f) {
	elem[4] += consts_1.twopi;
    }
    cea = rm * cv / elem[0] + elem[1];
/* Computing 2nd power */
    d__1 = elem[1];
    sea = rm * sv / (elem[0] * sqrt(1.f - d__1 * d__1));
    ea = atan2(sea, cea);
/* 	if(ea.lt.0.) ea = ea + twopi */
    per = elem[0] * sqrt(elem[0] / *emu);
    elem[5] = *t - (ea - elem[1] * sea) * per;
    per = consts_1.twopi * per;
L3:
    if (elem[5] < 0.) {
	elem[5] += per;
	goto L3;
    }
    ret_val = 0;
    return ret_val;
} /* ilmnts_ */


/* 	--------------------------------------------------- */
/* 	NK: changed to a function */
/* 	subroutine rrdot(elem,emu,t,r,rd) */
integer irrdot_(doublereal *elem, doublereal *emu, doublereal *t, doublereal *
	r__, doublereal *rd)
{
    /* System generated locals */
    integer ret_val;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal ea, em, en, rz[3];
    extern integer kepler_(doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int rotate_(doublereal *, doublereal *, integer *,
	     doublereal *);
    static doublereal edo;
    static integer kep;
    static doublereal rzd[3], srt;


/*       given r and rdot compute orbital elements */


/* 	calculation of eccentric anamoly from mean motion */
/* 	en = mean motion */

    en = sqrt(*emu / elem[0]) / elem[0];
    em = en * (*t - elem[5]);
/* 	NK: changed this as kepler is now a function returning int */
    kep = kepler_(&em, &elem[1], &ea);
    if (kep != 0) {
	ret_val = 1;
	return ret_val;
    }

/* 	calculate rectangular coordinates in orbital plane */

/* Computing 2nd power */
    d__1 = elem[1];
    srt = sqrt(1.f - d__1 * d__1);
    edo = en / (1.f - elem[1] * cos(ea));
    rz[0] = elem[0] * (cos(ea) - elem[1]);
    rz[1] = elem[0] * srt * sin(ea);
    rz[2] = 0.f;
    rzd[0] = -elem[0] * edo * sin(ea);
    rzd[1] = elem[0] * srt * edo * cos(ea);
    rzd[2] = 0.f;

/* 	now rotate to the ecliptic coordiante system */

    d__1 = -elem[4];
    rotate_(r__, &d__1, &c__3, rz);
    d__1 = -elem[4];
    rotate_(rd, &d__1, &c__3, rzd);
    d__1 = -elem[2];
    rotate_(rz, &d__1, &c__1, r__);
    d__1 = -elem[2];
    rotate_(rzd, &d__1, &c__1, rd);
    d__1 = -elem[3];
    rotate_(r__, &d__1, &c__3, rz);
    d__1 = -elem[3];
    rotate_(rd, &d__1, &c__3, rzd);
    ret_val = 0;
    return ret_val;
} /* irrdot_ */



/* 	--------------------------------------------------- */

/* 	NK: this was writing to stdout and messing up lc */
/* 	so I converted it to a function so that non-zero */
/* 	result flags an error */
/* 	subroutine kepler(em,e,ea) */
integer kepler_(doublereal *em, doublereal *e, doublereal *ea)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    double d_mod(doublereal *, doublereal *), sin(doublereal), cos(doublereal)
	    ;

    /* Local variables */
    static integer i__;
    static doublereal del, emr;



/* 	iterative solution of kepler's equation */

    emr = d_mod(em, &consts_1.twopi);
    *ea = emr;
    for (i__ = 1; i__ <= 99; ++i__) {
	del = (*ea - *e * sin(*ea) - emr) / (1.f - *e * cos(*ea));
	*ea -= del;
/* 		if(dabs(del).le.(1.d-7*(1. + dabs(ea)))) return */
	if (abs(del) <= (abs(*ea) + 1.f) * 1e-7) {
	    ret_val = 0;
	    return ret_val;
	}
/* L1: */
    }
/* 	write(*,3) em, e */
/* 	3	format(' Kepler failed for ',2f15.5) */
/* 	2	return */
    ret_val = 1;
    return ret_val;
} /* kepler_ */


/* 	--------------------------------------------------- */

/* Subroutine */ int rotate_(doublereal *r__, doublereal *a, integer *k, 
	doublereal *rz)
{
    /* Initialized data */

    static doublereal aold = 0.;
    static doublereal c__ = 1.;
    static doublereal s = 0.;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer i__, j;


/* 	rotate vector rz by angle a about the k axis */
/* 	k = 1, 2, 3 for x, y, z respectively */

    if (*a != aold) {
	c__ = cos(*a);
	s = sin(*a);
	aold = *a;
    }
    i__ = *k % 3 + 1;
    j = i__ % 3 + 1;
    r__[(i__1 = *k - 1) < 3 && 0 <= i__1 ? i__1 : s_rnge("r", i__1, "rotate_",
	     (ftnlen)161)] = rz[(i__2 = *k - 1) < 3 && 0 <= i__2 ? i__2 : 
	    s_rnge("rz", i__2, "rotate_", (ftnlen)161)];
    r__[(i__1 = i__ - 1) < 3 && 0 <= i__1 ? i__1 : s_rnge("r", i__1, "rotate_"
	    , (ftnlen)162)] = c__ * rz[(i__2 = i__ - 1) < 3 && 0 <= i__2 ? 
	    i__2 : s_rnge("rz", i__2, "rotate_", (ftnlen)162)] + s * rz[(i__3 
	    = j - 1) < 3 && 0 <= i__3 ? i__3 : s_rnge("rz", i__3, "rotate_", (
	    ftnlen)162)];
    r__[(i__1 = j - 1) < 3 && 0 <= i__1 ? i__1 : s_rnge("r", i__1, "rotate_", 
	    (ftnlen)163)] = -s * rz[(i__2 = i__ - 1) < 3 && 0 <= i__2 ? i__2 :
	     s_rnge("rz", i__2, "rotate_", (ftnlen)163)] + c__ * rz[(i__3 = j 
	    - 1) < 3 && 0 <= i__3 ? i__3 : s_rnge("rz", i__3, "rotate_", (
	    ftnlen)163)];
    return 0;
} /* rotate_ */


/* 	--------------------------------------------------- */

/* Subroutine */ int vcopy_(doublereal *v1, doublereal *v2)
{

/* 	set v2 = v1 */

    v2[0] = v1[0];
    v2[1] = v1[1];
    v2[2] = v1[2];
    return 0;
} /* vcopy_ */


/* 	--------------------------------------------------- */

/* Subroutine */ int vsum_(doublereal *x, doublereal *y, doublereal *z__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer i__;


/* 	z = x + y */

    for (i__ = 1; i__ <= 3; ++i__) {
	z__[(i__1 = i__ - 1) < 3 && 0 <= i__1 ? i__1 : s_rnge("z", i__1, 
		"vsum_", (ftnlen)190)] = x[(i__2 = i__ - 1) < 3 && 0 <= i__2 ?
		 i__2 : s_rnge("x", i__2, "vsum_", (ftnlen)190)] + y[(i__3 = 
		i__ - 1) < 3 && 0 <= i__3 ? i__3 : s_rnge("y", i__3, "vsum_", 
		(ftnlen)190)];
/* L1: */
    }
    return 0;
} /* vsum_ */


/* 	--------------------------------------------------- */

/* Subroutine */ int vsub_(doublereal *x, doublereal *y, doublereal *z__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer i__;


/* 	z = x - y */

    for (i__ = 1; i__ <= 3; ++i__) {
	z__[(i__1 = i__ - 1) < 3 && 0 <= i__1 ? i__1 : s_rnge("z", i__1, 
		"vsub_", (ftnlen)204)] = x[(i__2 = i__ - 1) < 3 && 0 <= i__2 ?
		 i__2 : s_rnge("x", i__2, "vsub_", (ftnlen)204)] - y[(i__3 = 
		i__ - 1) < 3 && 0 <= i__3 ? i__3 : s_rnge("y", i__3, "vsub_", 
		(ftnlen)204)];
/* L1: */
    }
    return 0;
} /* vsub_ */


/* 	--------------------------------------------------- */

/* Subroutine */ int vzero_(doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer i__;


/* 	x = 0 */

    for (i__ = 1; i__ <= 3; ++i__) {
	x[(i__1 = i__ - 1) < 3 && 0 <= i__1 ? i__1 : s_rnge("x", i__1, "vzer"
		"o_", (ftnlen)218)] = 0.f;
/* L1: */
    }
    return 0;
} /* vzero_ */


/* 	--------------------------------------------------- */

doublereal vdot_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    doublereal ret_val;


/* 	computer the dot project of 2 vectors x and y */

    ret_val = x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    return ret_val;
} /* vdot_ */


/* 	--------------------------------------------------- */

/* Subroutine */ int vxp_(doublereal *x, doublereal *y, doublereal *z__)
{

/* 	compute the vector cross product z = x cross y */

    z__[0] = x[1] * y[2] - x[2] * y[1];
    z__[1] = y[0] * x[2] - y[2] * x[0];
    z__[2] = x[0] * y[1] - x[1] * y[0];
    return 0;
} /* vxp_ */


/* 	--------------------------------------------------- */

doublereal vlen_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern doublereal vdot_(doublereal *, doublereal *);


/* 	return the length of vector x */

    ret_val = sqrt(vdot_(x, x));
    return ret_val;
} /* vlen_ */


/* 	--------------------------------------------------- */

doublereal stp_(doublereal *x, doublereal *y, doublereal *z__)
{
    /* System generated locals */
    doublereal ret_val;


/* 	compute the scaler triple product of x dotted into */
/* 	y cross z */

    ret_val = x[0] * (y[1] * z__[2] - y[2] * z__[1]) + x[1] * (z__[0] * y[2] 
	    - z__[2] * y[0]) + x[2] * (y[0] * z__[1] - y[1] * z__[0]);
    return ret_val;
} /* stp_ */


/* 	--------------------------------------------------- */


