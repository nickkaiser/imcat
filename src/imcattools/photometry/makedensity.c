#define	usage "\n\n\n\
NAME\n\
	makedensity --- bin catalogue into FITS image\n\
\n\
SYNOPSIS\n\
	makedensity r x1 x2 nx .....  [-v val] [-c]\n\
\n\
DESCRIPTION\n\
	'makedensity' reads a catalogue from stdin, which must\n\
	contain at least some N >=1 dimensional vector and sums the counts\n\
	of objects (or with the -v option sums some specified object value) into\n\
	bins in a floating point FITS image which is sent to stdout.\n\
\n\
	The first argument is the name of the coordinate vector.  This is then\n\
	followed by N triplets giving, starting with the fastest coordinate dimension,\n\
	the range of dimension to be mapped and the corresponding number of pixels in\n\
	the output image.\n\
\n\
	Use the -v option to sum some numerical object value named 'val'.  If val is a scalar\n\
	then the output image will have the same dimensionality as the coordinate\n\
	vector, but if val is a vector or matrix then the output image will be of higher\n\
	dimensionality and contain a set of images containing the the sums of the\n\
	various components of val. For example, if the input catalogue contains\n\
	a three dimensional coordinate r[3] = {x,y,z}, and a M1 x M2 matrix valued quantity\n\
	m[M2][M1], then the result is a 5-dimensional image f[N5][N4][N3][N2][N1]\n\
	with N5 = M2, N4 = M1 (and N3, N2, N1 given in the command line arguments).\n\
\n\
	With the -c option we assign the `charge' (v or unity) to four neighbouring\n\
	pixels.  This is done in such a way that if dx = (x2 - x1) / nx etc then a point\n\
	with x = ix + dx / dx is assigned entirely to the pixel with index ix.  The\n\
	model here is that the zeroth pixel extends from x = 0 to x = dx, etc.\n\
\n\
	The coordinate ranges are stored as FITS header records named\n\
		x0min, x0max, x1min, x1max, ....\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n\n"		


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../catlib/cat.h"
#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../utils/arrays.h"
#include "../../utils/args.h"
#include "../../utils/stats_stuff.h"


void		assigncharge(int level);
void		assigncharge_cic(int level, int xoff, double q);
void		setf(int level, int xoff, void *v);

static float	*f;
static int	*nx, fsize[MAX_DIMS + 1], vndim, *vdim, xdim, gxoff;
static double	*x, *x1, *x2, *dx, gq;
static void	*gv;
static char	*vname;

#define QSCALE_TINY	1.e-8

main(int argc, char *argv[])	
{
	char		*fitsfilename, *xname, *flag, limitstring[128];
	fitsheader	*fits;
	cathead		*inputcathead;
	object		*inputobject;
	item		*vitem, *xitem;
	int		vindex, xindex, dim, fndim, *fdim, npix, docic;
	int		i, xoff, ix;

	/* defaults */
	vname = NULL;
	docic = 0;
	
	/* start parsing args */
	argsinit(argc, argv, usage);
	if (nextargtype() != TEXT_ARG) {
		error_exit(usage);
	}
	xname = getargs();
	
	/* now read the catalogue header and see what kind of object 'xname' is */
	inputcathead = readcathead();			/* read the cat head */
	inputobject = newobject(inputcathead);		/* make the input object */
	connectobjecttocathead(inputobject);		/* obj addresses point back to cathead */
	allocobjectcontents(inputobject);		/* and allocate space for obj data */

	/* check the x-coordinate vector and find its dimension */
	xitem = getobjectitem(xname, inputcathead);
	if (xitem->itype != NUM_TYPE) {
		error_exit("makedensity: x must be numeric\n");
	}
	if (xitem->ndim != 1) {
		error_exit("makedensity: x must be 1-D vector\n");
	}
	xindex = getobjectitemindex(xname, inputobject);
	x = (double *) ((inputobject->addrlist)[xindex]);
	xdim = xitem->dim[0];

	/* now get the tripets 'x1 x2 nx' */
	x1 = (double *) calloc(xdim, sizeof(double));
	x2 = (double *) calloc(xdim, sizeof(double));
	nx = (int *) calloc(xdim, sizeof(int));
	dx = (double *) calloc(xdim, sizeof(double));
	for (dim = 0; dim < xdim; dim++) {
		x1[dim] = getargd();
		x2[dim] = getargd();
		nx[dim] = getargi();
		dx[dim] = (x2[dim] - x1[dim]) / nx[dim];
	}

	/* and loop over rest of the options */
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'v':
				vname = getargs();
				break;
			case 'c':
				docic = 1;
				break;
			default:
				argserror("unrecognized option", 0);
		}
	}

	/* check object val v if present */
	if (vname) {
		vitem = getobjectitem(vname, inputcathead);
		if (vitem->itype != NUM_TYPE) {
			error_exit("makedensity: v must be numeric\n");
		}
		vindex = getobjectitemindex(vname, inputobject);
		vndim = vitem->ndim;
		vdim = vitem->dim;
		gv = (inputobject->addrlist)[vindex];
	}
	
	/* contruct the fits header */
	fndim = xdim;
	if (vname) {
		if (vndim > 1 || vitem->dim[0] > 1) {
			fndim = xdim + vndim;
		}
	}
	fdim = (int *) calloc(fndim, sizeof(int));
	for (dim = 0; dim < xdim; dim++) {
		fdim[dim] = nx[dim];
	}
	if (fndim > xdim) {
		for (dim = 0; dim < vndim; dim++) {
			fdim[xdim + dim] = vitem->dim[dim];
		}
	}
	fits = newfitsheader(fndim, fdim, FLOAT_PIXTYPE);
	add_comment(argc, argv, fits);
	for (dim = 0; dim < xdim; dim++) {
		sprintf(limitstring, "x%dmin", dim);
		appendcomment(newnumericcomment(limitstring, x1[dim], NULL), fits);
		sprintf(limitstring, "x%dmax", dim);
		appendcomment(newnumericcomment(limitstring, x2[dim], NULL), fits);
	}
	writefitsheader(fits);

        /* figure out size of pixel, line, plane, cube etc of fits image in units of float */
        fsize[0] = 1;
        for (i = 1; i <= fits->ndim; i++) {
                fsize[i] = fsize[i-1] * fits->n[i-1];
        }

	/* allocate the data as a single block */
	f = (float *) calloc(fsize[fits->ndim], sizeof(float));
	if (!f) {
		error_exit("makedensity: failed to allocate memory for image\n");
	}


	/* loop over objects */
	while (readobject(inputobject)) {
		/* initialize pixel offset */
		gxoff = 0;
		gq = 1.0;
		if (docic) {
			assigncharge_cic(xdim - 1, 0, 1.0);
		} else {
			assigncharge(xdim - 1);
		}
	}

	/* write out the image */
	if (fits->opbyteorder != NATIVE_BYTE_ORDER) {
		byteswapline((void *) f, fsize[fits->ndim], pixsize(fits->extpixtype));
	}
	fwrite(f, sizeof(float), fsize[fits->ndim], stdout);
	exit(0);
}


void	assigncharge(int level)
{
	int	ix;

	ix = (int) floor((x[level] - x1[level]) / dx[level]);
	if ((ix >= 0) && (ix < nx[level])) {
		/* we are in allowed range */
		gxoff += ix * fsize[level];
		if (level) {
			assigncharge(level - 1);
		} else {
			if (!vname) {
				/* sum of counts */
				f[gxoff] += gq;
			} else {
				setf(vndim, gxoff, gv);
			}
		}
	}
}

void	assigncharge_cic(int level, int xoff, double q)
{
	int	ix, ix0;
	double	xx, qscale;

	/* scaled position of lower edge of cell */
	xx = (x[level] - x1[level]) / dx[level] - 0.5;
	/* lower pixel index */
	ix0 = (int) floor(xx);
	for (ix = ix0; ix <= ix0 + 1; ix++) {
		if ((ix >= 0) && (ix < nx[level])) {
			/* we are in allowed range */
			qscale = (ix == ix0 ? ix + 1 - xx : xx - ix0);
			if (qscale > QSCALE_TINY) {
				xoff += ix * fsize[level];
				if (level) {
					assigncharge_cic(level - 1, xoff, qscale * q);
				} else {
					gq = qscale * q;
					if (!vname) {
						/* sum of counts */
						f[xoff] += gq;
					} else {
						setf(vndim, xoff, gv);
					}
				}
				xoff -= ix * fsize[level];
			}
		}
	}
}

void	setf(int level, int xoff, void *v)
{
	int 		i, xshift;
	float		*fptr;

	level--;
	for (i = 0; i < vdim[vndim - level - 1]; i++) {
		xshift = i * fsize[level + xdim];
		if (level) {
			setf(level, xoff + xshift, ((void **) v)[i]);
		} else {
			f[xoff + xshift] += gq * ((double *) v)[i];
		}
	}
}

