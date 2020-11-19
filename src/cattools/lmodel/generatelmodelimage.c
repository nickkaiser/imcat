/*
 * generatelmodelimage.c
 */


#define usage "\n\
NAME\n\
	generatelmodelimage --- generate realisation of a lmodel as a FITS image\n\
\n\
SYNOPSIS\n\
	generatelmodelimage x1 x2 nx y1 y2 ny .... \n\
\n\
DESCRIPTION\n\
	'generatelmodelimage' reads from stdin a 'lmodel'\n\
	and computes a FITS image containing the realisation of the model\n\
	function on a grid of points spanning the rectangle bounded by\n\
	x = x1, x2; y = y1, y2; etc and with nx samples in x etc.\n\
	There must be one triplet of arguments for each dimension of x.\n\
	For example, for a 2-dimensional x, the command\n\
\n\
		generatelmodelimage 0 1 512 0 1 512\n\
\n\
	will generate an image of the model on the unit square with 512\n\
	pixels in each dimension.\n\
\n\
	If the dependent variable a is a scalar then the dimensionality\n\
	of the image is the same as that of the independent variable x.\n\
	For a 3-vector x[] for instance\n\
\n\
		f[iz][iy][ix] = sum_m a_m f_m(x)\n\
	with\n\
		x[0] = x1 + ix * (x2 - x1) / nx\n\
		x[1] = y1 + iy * (y2 - y1) / ny\n\
	etc.\n\
\n\
	If a is a matrix a[j0][j1]....[jn] then the image dimensionality\n\
	is the sum of the rank of a and the length of x, and the pixel\n\
	values are, for 2-vector x[] say\n\
\n\
		f[j0][j1]....[jn][iy][ix] = sum_m a_m[j0][j1]....[jn] f_m(x)\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "catlib/cat.h"
#include "imlib/fits.h"
#include "utils/error.h"
#include "utils/args.h"
#include "utils/lmodel.h"

void	*allocfitsarray(fitsheader *fits, int level);
int	writefitsarray(void *f, fitsheader *fits, int level);
int	addtofitsarray(void *fdst, void *fsrc, double fac, fitsheader *fits, int level);
int	generatefm(void *fm, int m, lmodel *themodel, fitsheader *fm_fits, int level);
int	aloop(item *theitem, void *a, void *f, void *fm, fitsheader *fm_fits, int level);



static	double	*x1, *x2, *dx, *x;


main(int argc, char *argv[])
{
	lmodel		*themodel;
	void		*a;
	int		m, i, *nx, fitsndim, *fitsdim, fm_ndim, *fm_dim;
	fitsheader	*fits, *fm_fits;
	void		*f, *fm;
	char		wcsstring[128];
	fitscomment	*thecomment;

	/* check argc */
	argsinit(argc, argv, usage);
	if (argc < 4) {
		error_exit(usage);
	}

	/* read the lmodel */
	themodel = readlmodel(stdin);

	/* get the image dimension arrays */
	x1 = (double *) calloc(themodel->xdim, sizeof(double));	
	x2 = (double *) calloc(themodel->xdim, sizeof(double));	
	dx = (double *) calloc(themodel->xdim, sizeof(double));	
	x = (double *) calloc(themodel->xdim, sizeof(double));	
	nx = (int *) calloc(themodel->xdim, sizeof(int));
	for (i = 0; i < themodel->xdim; i++) {
		x1[i] = getargd();
		x2[i] = getargd();
		nx[i] = getargi();
		dx[i] = (x2[i] - x1[i]) / nx[i];
	}	

	/* make the image dimension arrays */
	fitsndim = themodel->xdim;
	fitsndim += (themodel->aitem)->ndim;
	fitsdim = (int *) calloc(fitsndim, sizeof(int));	
	for (i = 0; i < themodel->xdim; i++) {
		fitsdim[i] = nx[i];
	}
	for (i = 0; i < (themodel->aitem)->ndim; i++) {
		fitsdim[fitsndim - i - 1] = (themodel->aitem)->dim[i];
	}
	fits = newfitsheader(fitsndim, fitsdim, FLOAT_PIXTYPE);
	fits->intpixtype = FLOAT_PIXTYPE;
	add_comment(argc, argv, fits);

	/* add the WCS header values */
	for (i = 0; i < themodel->xdim; i++) {
		sprintf(wcsstring, "CRVAL%d", i + 1);
		thecomment = newnumericcomment(wcsstring, x1[i], "");
		appendcomment(thecomment, fits);
		sprintf(wcsstring, "CRPIX%d", i + 1);
		thecomment = newnumericcomment(wcsstring, 0.0, "");
		appendcomment(thecomment, fits);
		sprintf(wcsstring, "CDELT%d", i + 1);
		thecomment = newnumericcomment(wcsstring, dx[i], "");
		appendcomment(thecomment, fits);
		sprintf(wcsstring, "CTYPE%d", i + 1);
		thecomment = newtextcomment(wcsstring, "PIXEL", "");
		appendcomment(thecomment, fits);
	}

	/* write the fits header */
	if (fits->n[fits->ndim-1] == 1) {
		fits->ndim--;	
		writefitsheader(fits);
		fits->ndim++;
	} else {
		writefitsheader(fits);
	}	

	/* generate the header for the fm[x] image */
	fm_ndim = themodel->xdim;
	fm_dim = nx;
	fm_fits = newfitsheader(fm_ndim, fm_dim, FLOAT_PIXTYPE);
	fm_fits->intpixtype = FLOAT_PIXTYPE;
		
	/* allocate the fm[x] image */
	fm = allocfitsarray(fm_fits, fm_fits->ndim);

	/* allocate image data f[I][x] */
	f = allocfitsarray(fits, fits->ndim);


	/* loop over modes m */
	for (m = 0; m < themodel->nmodes; m++) {
		generatefm(fm, m, themodel, fm_fits, fm_fits->ndim);
		aloop(themodel->aitem, themodel->a[m], f, fm, fm_fits, 0); 
	}

	writefitsarray(f, fits, fits->ndim);
	writefitstail(fits);
	exit(0);
}



/* recursive function to loop over elements of a[][]... */
int	aloop(item *theitem, void *a, void *f, void *fm, fitsheader *fm_fits, int level)
{
	int	i;

	if (level < theitem->ndim - 1) {
		for (i = 0; i < theitem->dim[level]; i++) {
			aloop(theitem, ((void **) a)[i], ((void **) f)[i], fm, fm_fits, level + 1);
		}
	} else {
		for (i = 0; i < theitem->dim[level]; i++) {
			addtofitsarray(((void **) f)[i], fm, ((double *) a)[i], fm_fits, fm_fits->ndim);
		}
	}
}


/* recursive function to allocate a fits array */
void	*allocfitsarray(fitsheader *fits, int level)
{
	void	*f;
	int	i;

	level--;
	if (level) {
		f = calloc(fits->n[level], sizeof(void *));
		for (i = 0; i < fits->n[level]; i++) {
			((void **)f)[i] = allocfitsarray(fits, level);
		}
	} else {
		f = calloc(fits->n[level], pixsize(fits->intpixtype));
	}
	return (f);
}


/* recursive function to write a fits array */
int	writefitsarray(void *f, fitsheader *fits, int level)
{
	int	i;

	level--;
	if (level) {
		for (i = 0; i < fits->n[level]; i++) {
			writefitsarray(((void **) f)[i], fits, level);
		}
	} else {
		writefitsline(f, fits);
	}
	return (1);
}


/* recursive function to generate f_m[x] */
int	generatefm(void *fm, int m, lmodel *themodel, fitsheader *fm_fits, int level)
{
	int	i;

	level--;
	for (i = 0; i < fm_fits->n[level]; i++) {
		x[level] = x1[level] + i * dx[level];
		if (level) {
			generatefm(((void **)fm)[i], m, themodel, fm_fits, level);
		} else {
			((float *) fm)[i] = lmodelfunc(themodel, m, x);
		}
	}
}

/* recursive function to add fac times fsrc to fdst */
int	addtofitsarray(void *fdst, void *fsrc, double fac, fitsheader *fits, int level)
{
	int	i;

	level--;
	if (level) {
		for (i = 0; i < (fits->n)[level]; i++) {
			addtofitsarray(((void **) fdst)[i], ((void **) fsrc)[i], fac, fits, level);
		}
	} else {
		switch (fits->intpixtype) {
			case DBL_PIXTYPE:
				for (i = 0; i < (fits->n)[level]; i++) {
					((double *) fdst)[i] += fac * ((double *) fsrc)[i];
				}
				break;
			case FLOAT_PIXTYPE:
				for (i = 0; i < (fits->n)[level]; i++) {
					((float *) fdst)[i] += fac * ((float *) fsrc)[i];
				}
				break;
			case INT_PIXTYPE:
				for (i = 0; i < (fits->n)[level]; i++) {
					((int *) fdst)[i] += fac * ((int *) fsrc)[i];
				}
				break;
			case SHORT_PIXTYPE:
				for (i = 0; i < (fits->n)[level]; i++) {
					((short *) fdst)[i] += fac * ((short *) fsrc)[i];
				}
				break;
			case UCHAR_PIXTYPE:
				for (i = 0; i < (fits->n)[level]; i++) {
					((unsigned char *) fdst)[i] += fac * ((unsigned char *) fsrc)[i];
				}
				break;
			default:
				error_exit("addtofitsarray : illegal pixel type\n");
		}
	}
}

