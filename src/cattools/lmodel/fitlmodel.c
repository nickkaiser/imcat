/*
 * fitlmodel.c
 */


#define usage "\n\
NAME\n\
	fitlmodel --- fit for linear superposition of mode function\n\
\n\
SYNOPSIS\n\
	fitlmodel xname aname modeldefinition [options...] \n\
\n\
DESCRIPTION\n\
	'fitlmodel' reads from stdin a catalogue containing at least a\n\
	'postion vector' x called xname and some other variable a called aname,\n\
	whch may be a scalar, vector or matrix of arbitrary rank,\n\
	and fits for a model of a(x) as a linear superposition\n\
	of a set of mode functions:\n\
\n\
		a(x) = sum_M a_M f_M(x)\n\
\n\
	where the mode function coefficients a_M have the same\n\
	dimensionality as a.\n\
\n\
	The 'modeldefinition' is a combination of arguments which can be\n\
	one of:\n\
\n\
	-p lmin lmax\n\
		Fit a polynomial. The mode functions are labelled by a set\n\
		of indices p[] with same length as x[], the functions are\n\
\n\
			f_p = x0^p0 x1^p1 .... = product x_i^p_i\n\
\n\
		and the order l = sum p_i lies in the inclusive interval lmin-lmax.\n\
		An alternative parameterisation of the indices is in terms of\n\
		the order array l[i] = l - sum_i=0^i-1 p[i], in terms of which the\n\
		p-indices are p[i] = l[i] - l[i+1].\n\
\n\
	-z nmin nmax\n\
		Fit for Zernike polynomials of order min through nmax\n\
		as defined in Born and Wolf.\n\
\n\
	-f kmin kmax lbox\n\
		Fit a sum of Fourier modes labelled by compound index m = {k[], i}\n\
\n\
			f_m = cos(2 PI k.x / lbox)		i = 0\n\
			f_m = sin(2 PI k.x / lbox)		i = 1\n\
\n\
		where the modulus of the integerised wave number k lies in the\n\
		range kmin - kmax inclusive.\n\
\n\
	-F kmax lbox\n\
		As for -f option save that all modes in hypercube with edges +- kmax.\n\
\n\
\n\
	Other options:\n\
\n\
	-c	# generate covariance matrix 'covar'\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "catlib/cat.h"
#include "utils/lu.h"
#include "utils/error.h"
#include "utils/args.h"
#include "utils/lmodel.h"



main(int argc, char *argv[])
{
	/* matrix stuff */
	int 	*indx, m, n, d;
	double	**A, **B, det;
	double	**aa, **aatranspose;

	/* catalogue stuff */
	cathead	*ipcat;
	object	*ipobj;
	item	*xitem, *aitem;
	int	ix, ia;

	/* model stuff */
	lmodel	*themodel;
	char	*xname, *aname, *modeltypeflag, *optflag;
	int	nmin, nmax, lmin, lmax, modeltype, kmin, kmax, makecovar;
	double	fm, fn, *x, lbox;
	void	*a;

	/* defaults */
	makecovar = 0;

	/* parse first 2 args */
	if (argc < 4) {
		error_exit(usage);
	}
	argsinit(argc, argv, usage);
	xname = getargs();
	aname = getargs();

	/* read cat head, create dummy oject and get indices*/
	ipcat = readcathead();
	ipobj = newobject(ipcat);
        connectobjecttocathead(ipobj);
        allocobjectcontents(ipobj);
	xitem = getobjectitem(xname, ipcat);
	aitem = getobjectitem(aname, ipcat);
	x = (double *) (xitem->addr);
	a = aitem->addr;

	modeltypeflag = getflag();
	switch (modeltypeflag[0]) {
		case 'p':
			modeltype = POLYNOMIAL_LMODEL;
			lmin = getargi();
			lmax = getargi();
			themodel = newpolylmodel(xitem, aitem, lmin, lmax);
			break;
		case 'z':
			modeltype = ZERNIKE_LMODEL;
			nmin = getargi();
			nmax = getargi();
			themodel = newzernikelmodel(xitem, aitem, nmin, nmax);
			break;
		case 'f':
			modeltype = FOURIER_LMODEL;
			kmin = getargi();
			kmax = getargi();
			lbox = getargd();
			themodel = newfourierlmodel(xitem, aitem, kmin, kmax, 0, lbox);
			break;
		case 'F':
			modeltype = FOURIER_LMODEL;
			kmax = getargi();
			lbox = getargd();
			themodel = newfourierlmodel(xitem, aitem, 0, kmax, 1, lbox);
			break;
		default:
			error_exit("fitlmodel: unrecognised model definition flag\n");
	}
	while (optflag = getflag()) {
		switch (optflag[0]) {
			case 'c':
				makecovar = 1;
				break;
			default:
				error_exit("fitlmodel: unrecognised option\n");
		}
	}

	/* copy previous comments and add history string */
	(themodel->cat)->commentbase = ipcat->commentbase;
	addargscomment(argc, argv, themodel->cat);

	/* allocate arrays for linear algebra */
	A = (double **) calloc(themodel->nmodes, sizeof(double *));
	for (m = 0; m < themodel->nmodes; m++) {
        	A[m] = (double *) calloc(themodel->nmodes, sizeof(double));
	}
        indx = (int *) calloc(themodel->nmodes, sizeof(int));

	/* accumulate arrays */
	while (readobject(ipobj)) {
		for (m = 0; m < themodel->nmodes; m++) {
			fm = lmodelfunc(themodel, m, x);
			addtomatrix((themodel->a)[m], a, fm, aitem->ndim, aitem->dim, 0);
			for (n = 0; n < themodel->nmodes; n++) {
				fn = lmodelfunc(themodel, n, x);
				A[m][n] += fm * fn;
			}
		}
	}


	myludcmp(A, themodel->nmodes, indx, &det);
	flatten_a(themodel);
	for (d = 0; d < themodel->asize; d++) { 
		mylubksb(A, themodel->nmodes, indx, (themodel->aflat)[d]);
	}
	unflatten_a(themodel);
	/* create the covariance matrix */
	if (makecovar) {
		themodel->hascovar = 1;
		themodel->covar = (double **) calloc(themodel->nmodes, sizeof(double *));
		for (m = 0; m < themodel->nmodes; m++) {
			(themodel->covar)[m] = (double *) calloc(themodel->nmodes, sizeof(double));
		}
		invertmatrix(A, themodel->covar, themodel->nmodes);
	}

	writelmodel(themodel, stdout);
	exit(0);
}


