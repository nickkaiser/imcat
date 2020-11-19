#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define REAL double
#include "triangle.h"
#include "utils/args.h"
#include "utils/iostream.h"
#include "utils/ipbuff.h"
#include "catlib/cat.h"
#include "imlib/fits.h"
#include "utils/error.h"
#include "utils/arrays.h"
#include "painttriangle.h"
#include "fixedges.h"
#include "average.h"

#define usage "\nNAME\n\
	triangulatecat - Delauney tesselation of a catalogue\n\
\n\
SYNOPSIS\n\
	triangulatecat srccat [-c | -f x1 x2 Nx y1 y2 Ny] [-m | -M]\n\
\n\
DESCRIPTION\n\
	Triangulatecat performs Delauney tesselation using Shewchuk's code.\n\
	Input catalog should contain a 2-vector x[2] and scalar or 1-D vector f[N].\n\
\n\
	Default is to output the triangles x[3][2], f[3][N].\n\
	With -c option we output an average cat x[2] = xbar[2], f[N] = fbar[N].\n\
	Ditto with -c -m options.\n\
	With -c -M options we output a median cat x[2] = xbar[2], f[N] = fmedian[N].\n\
	With -f ... option we output a Nx by Ny FITS image consisting of\n\
	triangular linear ramp segments, or, with -m or -M options the triangles\n\
	are uniform and painted with the mean or median f[] value respectively.\n\
\n\
	Uses iostream utilities so use '-' for standard input, 'somecommand |' to\n\
	read from a process etc.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"


/* output mode */
#define OP_TRIANGLES	0
#define OP_AVGCAT	1
#define OP_FITS		2

double	mean(double **y, int j);
double	median(double **y, int j);

main (int argc, char *argv[])
{
	char			*switches = "zQ";
	struct triangulateio 	*in, *out;
	int			v, nvert, t, ntriangles, c, ncorners;
	iostream		*ipstream;
	char			*catfilename, *xname, *fname, *flag;
	cathead			*cath;
	object			*obj, *objbase;
	item			*fitem, *xitem;
	int			findex, xindex, fdim, nobj, favgmode, opmode, Nx, Ny, i, ix, iy;
	void			*f, *x;
	double			x1, x2, y1, y2;
	double			*xt[3], *ft[3];
	fitsheader		*fits;
	float			***fim;
	int			fitsndim, fitsdim[MAX_FITS_DIM];

	/* defaults */
	xname = "x";
	fname = "f";
	opmode = OP_TRIANGLES;
	favgmode = F_INTERP;

	/* parse args */
	argsinit(argc, argv, usage);
	if (nextargtype() == FLAG_ARG) {
		error_exit(usage);
	}
	ipstream = openiostream(getargs(), "r");
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'c':
				opmode = OP_AVGCAT;
				favgmode = F_MEAN;
				break;
			case 'f':
				opmode = OP_FITS;
				favgmode = F_INTERP;
				x1 = getargd();
				x2 = getargd();
				Nx = getargi();
				y1 = getargd();
				y2 = getargd();
				Ny = getargi();
				break;
			case 'm':
				favgmode = F_MEAN;
				break;
			case 'M':
				favgmode = F_MEDIAN;
				break;
			default:
				error_exit(usage);
		}
	}
	if ((opmode == OP_TRIANGLES && favgmode) || (opmode == OP_AVGCAT && favgmode == F_INTERP)) {
		error_exit("triangulatecat: illegal combination of options\n");
	}

	/* now read the catalogue header */
	setcatipf(ipstream->f);
	cath = readcathead();			/* read the cat head */
	obj = newobject(cath);			/* make an input object */
	connectobjecttocathead(obj);		/* obj addresses point back to cathead */
	allocobjectcontents(obj);		/* and allocate space for obj data */

	/* check the x-coordinate vector and find its dimension */
	xitem = getobjectitem(xname, cath);
	if (xitem->itype != NUM_TYPE) {
		error_exit("triangulatecat: x must be numeric\n");
	}
	if (xitem->ndim != 1 || xitem->dim[0] != 2) {
		error_exit("triangulatecat: x must be 2-vector\n");
	}
	xindex = getobjectitemindex(xname, obj);
	x = ((obj->addrlist)[xindex]);

	/* check the x-coordinate vector and find its dimension */
	xitem = getobjectitem(xname, cath);
	if (xitem->itype != NUM_TYPE) {
		error_exit("triangulatecat: x must be numeric\n");
	}
	if (xitem->ndim != 1 || xitem->dim[0] != 2) {
		error_exit("triangulatecat: x must be 2-vector\n");
	}
	xindex = getobjectitemindex(xname, obj);
	x = (double *) ((obj->addrlist)[xindex]);

	/* check f-object */
	fitem = getobjectitem(fname, cath);
	if (fitem->itype != NUM_TYPE) {
		error_exit("makedensity: f must be numeric\n");
	}
	if (fitem->ndim != 1) {
		error_exit("triangulatecat: f must be scalar or 1-D vector\n");
	}
	findex = getobjectitemindex(fname, obj);
	fdim = fitem->dim[0];
	f = (obj->addrlist)[findex];

	/* read the objects into a null terminated linked list */
	objbase = NULL;
	nobj = 0;
        while (1) {
 		obj = newobject(cath);
		allocobjectcontents(obj);
                if (!readobject(obj))
                        break;
		nobj++;
                obj->next = objbase;
                objbase = obj;
        }
	
	/* create and fill the triangulateio objects */
	in  = (struct triangulateio *) calloc(1, sizeof(struct triangulateio));
	out = (struct triangulateio *) calloc(1, sizeof(struct triangulateio));
	nvert = nobj;
	in->numberofpoints = nvert;
	in->numberofpointattributes = fdim;
	in->pointlist = (REAL *) calloc(2 * in->numberofpoints, sizeof(REAL));
	in->pointattributelist = (REAL *) calloc(fdim * in->numberofpoints, sizeof(REAL));
	obj = objbase;
	for (v = 0; v < nvert; v++) {
		in->pointlist[2 * v] = ((double *) ((obj->addrlist)[xindex]))[0];
		in->pointlist[2 * v + 1] = ((double *) ((obj->addrlist)[xindex]))[1];
		for (i = 0; i < fdim; i++) {
			in->pointattributelist[fdim * v + i] = ((double *) ((obj->addrlist)[findex]))[i];
		}
		obj = obj->next;
	}

	/* do the biz */
	triangulate(switches, in, out, NULL);

	switch (opmode) {
		case OP_TRIANGLES:
			fitem->ndim = 2;
			fitem->dim[0] = 3;
			fitem->dim[1] = fdim;
			xitem->ndim = 2;
			xitem->dim[0] = 3;
			xitem->dim[1] = 2;
			obj = newobject(cath);			/* make an input object */
			allocobjectcontents(obj);		/* and allocate space for obj data */
			x = ((obj->addrlist)[xindex]);
			f = ((obj->addrlist)[findex]);
			addargscomment(argc, argv, cath);
			writecathead(cath);
			break;
		case OP_AVGCAT:
			obj = newobject(cath);			/* make an input object */
			allocobjectcontents(obj);		/* and allocate space for obj data */
			x = ((obj->addrlist)[xindex]);
			f = ((obj->addrlist)[findex]);
			addargscomment(argc, argv, cath);
			writecathead(cath);
			break;
		case OP_FITS:
			fitsdim[0] = Nx;
			fitsdim[1] = Ny;
			if (fdim > 1) {
				fitsndim = 3;
				fitsdim[2] = fdim;
			} else {
				fitsndim = 2;
			}
			fits = newfitsheader(fitsndim, fitsdim, FLOAT_PIXTYPE);
			add_comment(argc, argv, fits);
			appendcomment(newnumericcomment("x1", x1, NULL), fits);
			appendcomment(newnumericcomment("x2", x2, NULL), fits);
			appendcomment(newnumericcomment("y1", y1, NULL), fits);
			appendcomment(newnumericcomment("y2", y2, NULL), fits);
			writefitsheader(fits);
			fim = (float ***) calloc(fdim, sizeof(float **));
			for (i = 0; i < fdim; i++) {
				allocFloatArray(&(fim[i]), Nx, Ny);
				for (iy = 0; iy < Ny; iy++) {
					for (ix = 0; ix < Nx; ix++) {
						fim[i][iy][ix] = FLOAT_MAGIC;
					}
				}
			}
			break;
		default:
			exit(-1);
	}

			
	ntriangles = out->numberoftriangles;
	ncorners = out->numberofcorners;
	for (c = 0; c < 3; c++) {
		xt[c] = (double *) calloc(2, sizeof(double));
		ft[c] = (double *) calloc(fdim, sizeof(double));
	}
	for (t = 0; t < ntriangles; t++) {
		for (c = 0; c < 3; c++) {
			for (i = 0; i < 2; i++) {
				xt[c][i] = in->pointlist[2 * out->trianglelist[ncorners * t + c] + i];
			}
			for (i = 0; i < fdim; i++) {
				ft[c][i] = in->pointattributelist[fdim * out->trianglelist[ncorners * t + c] + i];
			}
		}
		switch (opmode) {
			case OP_TRIANGLES:
				for (c = 0; c < 3; c++) {
					for (i = 0; i < 2; i++) {
						((double **) x)[c][i] = xt[c][i];
					}
					for (i = 0; i < fdim; i++) {
						((double **) f)[c][i] = ft[c][i];
					}
				}
				writeobject(obj);
				break;
			case OP_AVGCAT:
				for (i = 0; i < 2; i++) {
					((double *) x)[i] = mean(xt, i);
				}
				for (i = 0; i < fdim; i++) {
					switch (favgmode) {
						case F_MEAN:
							((double *) f)[i] = mean(ft, i);
							break;
						case F_MEDIAN:
							((double *) f)[i] = median(ft, i);
							break;
						default:
							error_exit("triangulatecat: illegal average mode\n");
					}
				}
				writeobject(obj);
				break;
			case OP_FITS:
				painttriangle(fim, fdim, x1, x2, Nx, y1, y2, Ny, xt, ft, favgmode);
				break;
			default:
				exit(-1);
		}
	}

	if (opmode != OP_FITS) {
		closeiostream(ipstream);
		exit(0);
	}

	for (i = 0; i < fdim; i++) {
		fixedges(fim[i], Nx, Ny);
		fixholes(fim[i], Nx, Ny);
		writefitsplane((void **) (fim[i]), fits);
	}
	writefitstail(fits);
	exit(0);
}



