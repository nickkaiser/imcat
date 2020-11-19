/* 
 * lmodel.c - functions to create and manipulate linear fit models ('lmodels')
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "catlib/cat.h"
#include "utils/error.h"
#include "lmodel.h"

/* create a modeltype independent lmodel structure */
lmodel *newlmodel(item *xitem, item *aitem)
{
	lmodel	*themodel;
	int	d;

	/* check the type of xitem, atime */
	if (xitem->itype != NUM_TYPE || xitem->ndim != 1 || aitem->itype != NUM_TYPE) {
		error_exit("fitlmodel: non-numerical items or non-vector position\n");
	}

	/* allocate the lmodel */
	themodel = (lmodel *) calloc(1, sizeof(lmodel));

	/* create new items */
	themodel->xitem = (item *) calloc(1, sizeof(item));
	themodel->aitem = (item *) calloc(1, sizeof(item));
	
	/* copy over contents of xitem, aitem */
	*(themodel->xitem) = *xitem;
	*(themodel->aitem) = *aitem;

	/* make copies of some useful stuff */
	themodel->xdim = xitem->dim[0];
	themodel->xname = xitem->name;
	themodel->aname = aitem->name;

	/* default origin */
	themodel->xorigin = (double *) calloc(themodel->xdim, sizeof(double));
	for (d = 0; d < themodel->xdim; d++) {
		(themodel->xorigin)[d] = 0.0;
	}

	/* size of a in doubles */
	themodel->asize = 1;
	for (d = 0; d < (themodel->aitem)->ndim; d++) {
		themodel->asize *= (aitem->dim)[d];
	}
	
	themodel->cat = (cathead *) calloc(1, sizeof(cathead));
	return(themodel);
}

/* create a full polymodel l-model structure */
lmodel *newpolylmodel(item *xitem, item *aitem, int lmin, int lmax)
{
	lmodel	*themodel;
	int	d, l, m;

	/* create generic lmodel */
	themodel = newlmodel(xitem, aitem);

	/* set the modeltype */
	themodel->modeltype = POLYNOMIAL_LMODEL;

	/* allocate space for index vector used by polyloop */
	themodel->pp = (int *) calloc(themodel->xdim, sizeof(int));

	/* set the order limits */
	themodel->lmin = lmin;
	themodel->lmax = lmax;

	/* compute nmodes */
	themodel->nmodes = 0;
	for (l = themodel->lmin; l <= themodel->lmax; l++) {
		polyloop(themodel, l, 0);
	}

	/* allocate the p-array, l-array */
	themodel->p = (double **) calloc(themodel->nmodes, sizeof(double *));
	themodel->l = (double **) calloc(themodel->nmodes, sizeof(double *));
	for (m = 0; m < themodel->nmodes; m++) {
		(themodel->p)[m] = (double *) calloc(themodel->xdim, sizeof(double));
		(themodel->l)[m] = (double *) calloc(themodel->xdim, sizeof(double));
	}

	/* allocate a-array */
	themodel->a = (void **) calloc(themodel->nmodes, sizeof(void *));
	for (m = 0; m < themodel->nmodes; m++) {
		allocitemcontents(themodel->aitem, themodel->a + m, 0);
	}

	/* assign p, l arrays */
	themodel->nmodes = 0;
	for (l = themodel->lmin; l <= themodel->lmax; l++) {
		polyloop(themodel, l, 0);
	}
	return(themodel);
}



/* recursive function to loop over all modes of order l */
int	polyloop(lmodel *themodel, int l, int i)
{
	int 	p, pmin, pmax, psum, j;

	psum = 0;
	for (j = 0; j < i; j++) {
		psum += (themodel->pp)[j];
	}
	if (i < themodel->xdim - 1) {
		for (p = 0; p <= l - psum; p++) {
			(themodel->pp)[i] = p;
			polyloop(themodel, l, i+1);
		} 
	} else {
		(themodel->pp)[i] = l - psum;	
		if (themodel->p) {
			psum = 0;
			for (j = 0; j < themodel->xdim; j++) {
				(themodel->l)[themodel->nmodes][j] = (double) (l - psum);
				(themodel->p)[themodel->nmodes][j] = (double) ((themodel->pp)[j]);
				psum += (themodel->pp)[j];
			}
		}
		(themodel->nmodes)++;
	}
	return(1);
}


/* create a full zernike l-model structure */
lmodel *newzernikelmodel(item *xitem, item *aitem, int nmin, int nmax)
{
	lmodel	*themodel;
	int	m, n, mm;

	/* check that x is a 2 vector */
	if (xitem->ndim != 1 || xitem->dim[0] != 2) {
		error_exit("newzernikelmodel : x must be a two vector\n");
	}

	/* create generic lmodel */
	themodel = newlmodel(xitem, aitem);
	themodel->modeltype = ZERNIKE_LMODEL;

	/* set the order limits */
	themodel->nmin = nmin;
	themodel->nmax = nmax;

	/* compute nmodes */
	themodel->nmodes = 0;
	for (n = nmin; n <= nmax; n++) {
		themodel->nmodes += n + 1;
	}

	/* allocate the n, m arrays */
	themodel->n = (double *) calloc(themodel->nmodes, sizeof(double));
	themodel->m = (double *) calloc(themodel->nmodes, sizeof(double));

	/* allocate a-array */
	themodel->a = (void **) calloc(themodel->nmodes, sizeof(void *));
	for (m = 0; m < themodel->nmodes; m++) {
		allocitemcontents(themodel->aitem, themodel->a + m, 0);
	}

	/* assign n, m arrays */
	mm = 0;
	for (n = themodel->nmin; n <= themodel->nmax; n++) {
		for (m = -n; m <= n; m += 2) {
			themodel->n[mm] = (double) n;
			themodel->m[mm] = (double) m;
			mm++;
		}
	}

	/* generate R array */
	makezernikeR(themodel);

	return(themodel);
}


/* generate the R array */
int	makezernikeR(lmodel *themodel)
{
	int	mm, n, m, mabs, s;

	themodel->R = (double **) calloc(themodel->nmodes, sizeof(double *));
	for (mm = 0; mm < themodel->nmodes; mm++) {
		n = (int) themodel->n[mm];
		m = (int) themodel->m[mm];
		mabs = abs(m);
		themodel->R[mm] = (double *) calloc(1 + (n - mabs) / 2, sizeof(double));
		for (s = 0; s <= (n - mabs) / 2; s++) {
			themodel->R[mm][s] = bang(n - s) / (bang(s) * bang((n + mabs) / 2 - s) * bang((n - mabs) / 2 - s));
			if (2 * (s / 2) != s) {
				themodel->R[mm][s] *= -1.0;
			}
		}
	}
}


/* create a full fourier l-model structure */
lmodel *newfourierlmodel(item *xitem, item *aitem, int kmin, int kmax, int boxlimits, double lbox)
{
	lmodel	*themodel;
	int 	m;

	/* create generic lmodel */
	themodel = newlmodel(xitem, aitem);
	/* set the modeltype */
	themodel->modeltype = FOURIER_LMODEL;

	/* allocate space for index vector used by fourierloop */
	themodel->kk = (int *) calloc(themodel->xdim, sizeof(int));

	/* set the order limits */
	themodel->kmin = kmin;
	themodel->kmax = kmax;
	themodel->boxlimits = boxlimits;
	themodel->lbox = lbox;

	/* compute nmodes */
	themodel->nmodes = 0;
	fourierloop(themodel, 0);

	/* allocate the k-array, i-array */
	themodel->k = (double **) calloc(themodel->nmodes, sizeof(double *));
	themodel->i = (double *) calloc(themodel->nmodes, sizeof(double));
	for (m = 0; m < themodel->nmodes; m++) {
		(themodel->k)[m] = (double *) calloc(themodel->xdim, sizeof(double));
	}

	/* allocate a-array */
	themodel->a = (void **) calloc(themodel->nmodes, sizeof(void *));
	for (m = 0; m < themodel->nmodes; m++) {
		allocitemcontents(themodel->aitem, themodel->a + m, 0);
	}

	/* assign k, i arrays */
	themodel->nmodes = 0;
	fourierloop(themodel, 0);

	return(themodel);
}


/* recursive function to loop over all fourier modes */
int	fourierloop(lmodel *themodel, int i)
{
	int 	j, isfirst, isnotorigin, ii;
	double	kk;


	/* isfirst is set if the previous components pp[j<i] are all zero */
	isfirst = 1;
	for (j = 0; j < i; j++) {
		if (themodel->kk[j]) {
			isfirst = 0;
			break;
		}
	}

	/* we loop over standard set of modes where first non-zero component of k must be non-negative */
	for (themodel->kk[i] = (isfirst ? 0 : -themodel->kmax); themodel->kk[i] <= themodel->kmax; (themodel->kk[i])++) {
		if (i < themodel->xdim - 1) {
			fourierloop(themodel, i + 1);
		} else {
			kk = 0.0;
			for (j = 0; j < themodel->xdim; j++) {
				kk += themodel->kk[j] * themodel->kk[j];
			}
			isnotorigin = (themodel->kk[j-1] == 0 && isfirst ? 0 : 1);
			if ((kk >= themodel->kmin * themodel->kmin) && (themodel->boxlimits || (kk <= themodel->kmax * themodel->kmax))) {
				/* if the k, i array are allocated we assign them */
				for (ii = 0; ii <= isnotorigin; ii++) {
					if (themodel->k) {
						for (j = 0; j < themodel->xdim; j++) {
							themodel->k[themodel->nmodes][j] = (double) (themodel->kk[j]);
						}
						themodel->i[themodel->nmodes] = (double) ii;
					}				
					themodel->nmodes++;
				}
			}
		}
	}
	return(1);
}




/* the mode functions */
double	lmodelfunc(lmodel *themodel, int m, double *x)
{
	/* general stuff */
	double	result;
	int	d;
	/* zernike stuff */	
	double	phi, r, R;
	int	in, im, mabs, s;
	/* fourier stuff */
	double	kx;

	switch (themodel->modeltype) {
		case POLYNOMIAL_LMODEL:
			result = 1.0;
			for (d = 0; d < themodel->xdim; d++) {
				result *= pow(x[d], (themodel->p)[m][d]);
			}
			return(result);
		case ZERNIKE_LMODEL:
			r = sqrt(x[0] * x[0] + x[1] * x[1]);
			if (r > 1.0) {
				return(0.0);
			}
			phi = atan2(x[1], x[0]);
			result = 0;
			in = (int) themodel->n[m];
			im = (int) themodel->m[m];
			mabs = abs(im);
			for (s = 0; s <= (in - mabs) / 2; s++) {
				result += themodel->R[m][s] * pow(r, in - 2 * s);
			}
			if (im < 0) {
				result *= sin(mabs * phi);
			} else {
				result *= cos(mabs * phi);
			}
			return(result);
		case FOURIER_LMODEL:
			kx = 0.0;
			for (d = 0; d < themodel->xdim; d++) {
				kx += themodel->k[m][d] * x[d];
			}
			if (themodel->i[m] > 0.0) {
				result = sin(2 * M_PI * kx / themodel->lbox);
			} else {
				result = cos(2 * M_PI * kx / themodel->lbox);
			}
			return(result);
		default:
			error_exit("lmodelfunc: unknown model type\n");
	}
}

/* functions to convert from a[m][][]... to aflat[j][m] and vice versa */
int	flatten_a(lmodel *themodel)
{
	int	j, m;

	if (!(themodel->aflat)) {
		themodel->aflat = (double **) calloc(themodel->asize, sizeof(double *));
		for (j = 0; j < themodel->asize; j++) {
			(themodel->aflat)[j] = (double *) calloc(themodel->nmodes, sizeof(double));
		}
	}
	for (m = 0; m < themodel->nmodes; m++) {
		themodel->j = 0;
		rflatten_a(themodel, (themodel->a)[m], m, 0);
	}
}


int	rflatten_a(lmodel *themodel, void *a, int m, int level)
{
	int	i;

	if (level < (themodel->aitem)->ndim - 1) {
		for (i = 0; i < ((themodel->aitem)->dim)[level]; i++) {
			rflatten_a(themodel, ((void **) a)[i], m, level + 1);
		}
	} else {
		for (i = 0; i < ((themodel->aitem)->dim)[level]; i++) {
			(themodel->aflat)[themodel->j++][m] = ((double *) a)[i];
		}
	}
}

int	unflatten_a(lmodel *themodel)
{
	int	m;

	for (m = 0; m < themodel->nmodes; m++) {
		themodel->j = 0;
		runflatten_a(themodel, (themodel->a)[m], m, 0);
	}
}


int	runflatten_a(lmodel *themodel, void *a, int m, int level)
{
	int	i;

	if (level < (themodel->aitem)->ndim - 1) {
		for (i = 0; i < ((themodel->aitem)->dim)[level]; i++) {
			runflatten_a(themodel, ((void **) a)[i], m, level + 1);
		}
	} else {
		for (i = 0; i < ((themodel->aitem)->dim)[level]; i++) {
			((double *) a)[i] = (themodel->aflat)[themodel->j++][m];
		}
	}
}


/* recursive function to add fac times asrc to a */
int	addtomatrix(void *adst, void *asrc, double fac, int ndim, int *dim, int level)
{
	int	i;

	if (level < ndim - 1) {
		for (i = 0; i < dim[level]; i++) {
			addtomatrix(((void **) adst)[i], ((void **) asrc)[i], fac, ndim, dim, level + 1);
		}
	} else {
		for (i = 0; i < dim[level]; i++) {
			((double *) adst)[i] += fac * ((double *) asrc)[i];
		}
	}
}


/* recursive function to zero a matrix */
int	zeromatrix(void *a, int ndim, int *dim, int level)
{
	int	i;

	if (level < ndim - 1) {
		for (i = 0; i < dim[level]; i++) {
			zeromatrix(((void **) a)[i], ndim, dim, level + 1);
		}
	} else {
		for (i = 0; i < dim[level]; i++) {
			((double *) a)[i] = 0.0;
		}
	}
}


/* factorial */
int	bang(int n)
{
	int result = 1;

	while (n > 1) {
		result *= n;
		n--;
	}
	return(result);
}

