#define	usage "\n\n\n\
NAME\n\
	getfitsval --- add values from a fits file to a catalogue\n\
\n\
SYNOPSIS\n\
	getfitsval fitsfile valname [options...]\n\
		-x xname	# spatial coordinate ('x')\n\
		-m magicval	# magic value (-32768)\n\
\n\
DESCRIPTION\n\
	'getfitsval' reads a catalogue from stdin, which must\n\
	contain at least some coordinate 'x' and adds an entry\n\
	named 'valname' with value derived from the FITS file\n\
	'fitsfile' and writes result to stdout.\n\
\n\
	The dimensionality of the new entry depends on the dimensionality\n\
	of the x-coordinate and the image f.  If these match then the new\n\
	entry is a scalar; for a two dimensional coordinate and a two\n\
	dimensional image for example, the output value is v = f[iy][ix] where\n\
		iy = (int) floor(x[1])\n\
		ix = (int) floor(x[0])\n\
	If the FITS image is of higher dimension, then the new value will\n\
	be a vector or matrix.  For example, with a 2-vector x and 3D image\n\
	f[iz][iy][ix] with dimensions Nz, Ny, Nx, the output value v will be a\n\
	vector of size v[Nz], with values\n\
		v[iz] = f[iz][iy][ix]\n\
	Similarly, if f has five dimensions say, and x is a 3-vector, then\n\
	v is a matrix of dimensions v[N5][N4].\n\
\n\
	By default, getfitsval looks for a spatial coordinate named 'x' but\n\
	you can substitute another name with the '-x' option.\n\
\n\
	Points which lie in pixels with MAGIC value or which\n\
	fall outside the image are assigned the value 'magicval'.\n\
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

void	setv(int level, int xoff, void *v);
static double magicval;
static	int	vndim, vdim[MAX_DIMS], fsize[MAX_DIMS + 1], xdim;
static	char	*f;


main(int argc, char *argv[])	
{
	char		*fitsfilename, *xname, defaultxname[2] = "x", *vname, *flag;
	fitsheader	*fits;
	FILE		*fitsf;
	cathead		*inputcathead, *outputcathead;		/* catalogue stuff... */
	object		*inputobject, *outputobject;
	item		*vitem, *xitem;
	int		vindex, idim, i, ix, xoff, nlines, cattype;
	void		*v;
	double		*x;


	/* defaults */
	xname = defaultxname;
	magicval = SHORT_MAGIC;
	
	/* parse args */
	argsinit(argc, argv, usage);
	if (argc < 3)
		error_exit(usage);
	fitsfilename = getargs();
	vname = getargs();
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'x':
				xname = getargs();
				break;
			case 'm':
				magicval = (double) getargf();
				break;
			default:
				error_exit(usage);
		}
	}

	/* read the fits file header */
	fitsf = fopen(fitsfilename, "r");
	if (!fitsf) {
		error_exit("getfitsval: unable to open fits file\n");
	}
	fits = readfitsheader(fitsf);

	/* set up inputcathead and inputobject */
	inputcathead = readcathead();			/* read the cat head */
	inputobject = newobject(inputcathead);		/* make the input object */
	connectobjecttocathead(inputobject);		/* obj addresses point back to cathead */
	allocobjectcontents(inputobject);		/* and allocate space for obj data */

	/* check the x-coordinate vector and find its dimension */
	xitem = getobjectitem(xname, inputcathead);
	if (xitem->itype != NUM_TYPE) {
		error_exit("getfitsval: x must be numeric\n");
	}
	if (xitem->ndim != 1) {
		error_exit("getfitsval: x must be 1-D vector\n");
	}
	xdim = xitem->dim[0];

	/* figure out the dimensions of v */
	if (xdim > fits->ndim) {
		error_exit("getfitsval: dimensionality of x exceeds that of image\n");
	}
	if (xdim == fits->ndim) {	/* v is a scalar */
		vndim = 1;
		vdim[0] = 1;
	} else {			/* v is a vector or matrix */
		vndim = fits->ndim - xdim;
		if (vndim >= MAX_DIMS) {
			error_exit("getfitsval: dimensions of created value too large\n");
		}
		for (idim = 0; idim < vndim; idim++) {
			vdim[idim] = fits->n[fits->ndim - idim - 1];
		}
	}

	/* construct the output cathead */
	outputcathead = (cathead *) calloc(1, sizeof(cathead));			/* new cathead */
	copyheaderinfo(outputcathead, inputcathead);				/* copy header stuff */
	addargscomment(argc, argv, outputcathead);		/* add history */
	copycontentinfo(outputcathead, inputcathead);		/* copy over pre-exisiting object items */
	getcatipfiletype(&cattype);
	setcatopfiletype(cattype);	

	/* construct the new item */
	vitem = newitembydimarray(vname, NUM_TYPE, vndim, vdim);
	addobjectitem(vitem, outputcathead);
	writecathead(outputcathead);				/* and write cathead out */			

	outputobject = newobject(outputcathead);		/* make the output object */
	vindex = getobjectitemindex(vname, outputobject);	/* get indices for new items */
	inheritcontents(outputobject, inputobject);		/* pre-existing output object addresses point back to inputobject */
	allocitemcontents(vitem, &((outputobject->addrlist)[vindex]), 0);	/* allocate space for new data */
	v = (outputobject->addrlist)[vindex];			/* get local handle for the new item. */

	/* now we get the handles to the input object items we will need */
	x = (double *) ((inputobject->addrlist)[getobjectitemindex(xname, inputobject)]);

	/* figure out nlines */
	nlines = 1;
	for (i = 1; i < fits->ndim; i++) {
		nlines *= fits->n[i];
	}

	/* figure out size of pixel, line, plane, cube etc of fits image */
	fsize[0] = pixsize(fits->intpixtype);
	for (i = 1; i <= fits->ndim; i++) {
		fsize[i] = fsize[i-1] * fits->n[i-1];
	}

	/* allocate the data as 1-D array of chars */
	f = (char *) calloc(fsize[fits->ndim], sizeof(char));
	/* and read the image data */
	for (i = 0; i < nlines; i++) {
		readfitsline((float *) (f + i * fsize[1]), fits);
	}

	while (readobject(inputobject)) {
		/* work out the offset in chars */
		xoff = 0;
		for (i = 0; i < xdim; i++) {
			ix = (int) floor(x[i]);
			if ((ix < 0) || (ix >= fits->n[i])) {	/* out of range */
				xoff = -1;
				break;
			} else {
				xoff += ix * fsize[i];
			}
		}
		setv(vndim, xoff, v);
		writeobject(outputobject);
	}
	exit(0);
}


void	setv(int level, int xoff, void *v)
{
	int 		i, xshift;
	float		*fptr;

	level--;
	for (i = 0; i < vdim[vndim - level - 1]; i++) {
		xshift = i * fsize[level + xdim];
		if (level) {
			if (xoff >= 0) {
				setv(level, xoff + xshift, ((void **) v)[i]);
			} else {
				setv(level, xoff, ((void **) v)[i]);
			}
		} else {
			if (xoff >= 0) {
				fptr = (float *) (f + xoff + xshift);
				if (*fptr == FLOAT_MAGIC) {
					((double *) v)[i] = magicval;	
				} else {
					((double *) v)[i] = (double) (*fptr);
				}
			} else {
				((double *) v)[i] = magicval;	
			}
		}
	}
}

