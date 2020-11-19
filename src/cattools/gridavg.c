/*
 * pastecats.c
 */

#define usage "\n\n\
NAME\n\
	gridavg - average objects binned onto a spatial grid\n\
\n\
SYNOPSIS\n\
	gridavg d [options....]\n\
		-x xname	# name for spatial coord ('x')\n\
		-r x1 x2 y1 y2	# range of spatial coords\n\
		-m		# take median\n\
\n\
DESCRIPTION\n\
	'gridavg' first reads a catalogue of objects from stdin\n\
	and then assigns them to cells in a grid with spacing d.\n\
	We then average the object values for each grid cell and\n\
	output the resulting catalogue to stdout.  The output catalogue\n\
	contains a leading column 'ncell' containing the count\n\
	of objects on the cell and following columns contain average\n\
	values for the numeric items in the input catalogue.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n"


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "../utils/error.h"
#include "../utils/fmedian.h"
#include "../catlib/cat.h"


typedef struct numobj {
	double 		*f;
	struct numobj 	*next;
} numobj;

#define	MAX_NUMBERS 10000

static int		nnum, domedian;
static double		*gaddr[MAX_NUMBERS], *fout;
static float		**fmed;

void    getaddresses(item *theitem, void *addr, int level);
void	average(int count, numobj *baseobj);


main(int argc, char *argv[])
{
	int		arg;
	cathead		*ipcat, *opcat;
	object		*ipobj;
	numobj		*thenumobj, *numobjbase = NULL, ***grid, *next;
	item		*ipitem, *opitem;
	char		*xname, defaultxname[2] = "x";
	int		xind, yind, needxrange, first, inum;
	double		x, y, d, x1, x2, y1, y2;
	int		ix, iy, nx, ny, count;

	/* defaults */
	xname = defaultxname;
	needxrange = 1;
	domedian = 0;

	if (argc < 2)
		error_exit(usage);
	if (1 != sscanf(argv[1], "%lf", &d))
		error_exit(usage);

	/* parse args */
	arg = 2;
	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch (argv[arg++][1]) {
			case 'x':
				xname = argv[arg++];
				break;
			case 'r':
				sscanf(argv[arg++], "%lf", &x1);
				sscanf(argv[arg++], "%lf", &x2);
				sscanf(argv[arg++], "%lf", &y1);
				sscanf(argv[arg++], "%lf", &y2);
				needxrange = 0;
				break;
			case 'm':
				domedian = 1;
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	/* read the cat header and create input object*/
	ipcat  = readcathead();
	ipobj = newobject(ipcat);

	/* create the output cathead */
	opcat = (cathead *) calloc(1, sizeof(cathead));
	copyheaderinfo(opcat, ipcat);
	addargscomment(argc, argv, opcat);

	/* create the output object items */
	opitem = newitem("ncell", NUM_TYPE, 1, 1);
	addobjectitem(opitem, opcat);
	ipitem = ipcat->objectitembase;
	nnum = 0;
	yind = 0;
	while (ipitem) {
		allocitemcontents(ipitem, &(ipitem->addr), 0);
		if (ipitem->itype == NUM_TYPE) {
			opitem = copyitem(ipitem);
			opitem->addr = ipitem->addr;
			addobjectitem(opitem, opcat);
			if (!strcmp(ipitem->name, xname)) {
				if ((ipitem->ndim != 1) || ((ipitem->dim)[0] != 2))
					error_exit("gridavg: spatial coord is not a 2-vector!\n");
				xind = nnum;
				yind = xind + 1;
			}
			getaddresses(ipitem, ipitem->addr, 0);
		}
		ipitem = ipitem->next;
	}
	if (nnum > MAX_NUMBERS)
		error_exit("gridavg: too many mumeric items in catalogue\n");
	connectobjecttocathead(ipobj);
	if (!yind)
		error_exit("gridavg: could not find spatial coord!\n");

	/* write output cathead */
	setcatopfiletype(BINARY_FILE_TYPE);	
	writecathead(opcat);

	/* read the objects into a linked list and get xrange if necessary */
	first = 1;
	while (readobject(ipobj)) {
		thenumobj = (numobj *) calloc(1, sizeof(numobj));
		thenumobj->f = (double *) calloc(nnum, sizeof(double));
		for (inum = 0; inum < nnum; inum++) {
			(thenumobj->f)[inum] = *(gaddr[inum]);
		}
		thenumobj->next = numobjbase;
		numobjbase = thenumobj;
		if (needxrange) {
			x = *(gaddr[xind]);
			y = *(gaddr[yind]);
			if (first) {
				x1 = x2 = x;
				y1 = y2 = y;
				first = 0;
			} else {
				x1 = (x < x1 ? x : x1);
				x2 = (x > x2 ? x : x2);
				y1 = (y < y1 ? y : y1);
				y2 = (y > y2 ? y : y2);
			}
		}
	}

	/* figure size of grid and allocate it */
	if (!(x2 > x1 && y2 > y1))
		error_exit("gridavg: bad coordinate range!\n");
	nx = (int) ceil((x2 - x1) / d);
	ny = (int) ceil((y2 - y1) / d);
	grid = (numobj ***) calloc(nx, sizeof(numobj **));
	for (ix = 0; ix < nx; ix++) {
		grid[ix] = (numobj **) calloc(ny, sizeof(numobj *));
	}

	/* install the objects ob the grid */
	thenumobj = numobjbase;
	while (thenumobj) {
		x = (thenumobj->f)[xind];
		y = (thenumobj->f)[yind];
		ix = (int) floor((x - x1) / d);
		iy = (int) floor((y - y1) / d);
		next = thenumobj->next;
		if (ix >= 0 && ix < nx && iy >= 0 && iy < ny) {
			thenumobj->next = grid[ix][iy];
			grid[ix][iy] = thenumobj;
		}
		thenumobj = next;
	}

	/* loop over grid cells, count objects and average */
	fout = (double *) calloc(nnum + 1, sizeof(double));
	if (domedian) {
		fmed = (float **) calloc(nnum, sizeof(float*));
	}
	for (ix = 0; ix < nx; ix++) {
		for (iy = 0; iy < ny; iy++) {
			thenumobj = grid[ix][iy];
			count = 0;
			while (thenumobj) {
				count++;		
				thenumobj = thenumobj->next;
			}
			if (count) {
				average(count, grid[ix][iy]);
				fwrite(fout, sizeof(double), nnum + 1, stdout);
			}
		}
	}
	exit(0);
}




void    getaddresses(item *theitem, void *addr, int level)
{
        int     i;

        if (level < theitem->ndim - 1) {
                for (i = 0; i < (theitem->dim)[level]; i++) {
                        getaddresses(theitem, *((void **) addr + i), level + 1);
                }
        } else {
                for (i = 0; i < (theitem->dim)[level]; i++) {
                        gaddr[nnum++] = (double *) addr + i;
                }
        }
}


void	average(int count, numobj *baseobj)
{
	int	inum, iobj;
	numobj	*theobj;

	fout[0] = (double) count;
	for (inum = 0; inum < nnum; inum++) {
		if (domedian) {
			fmed[inum] = (float *) calloc(count, sizeof(float));
		} else {
			fout[inum + 1] = 0.0;
		}
	}
	theobj = baseobj;
	iobj = 0;
	while (theobj) {
		for (inum = 0; inum < nnum; inum++) {
			if (domedian) {
				fmed[inum][iobj] = (float) ((theobj->f)[inum]);
			} else {
				fout[inum + 1] += (theobj->f)[inum];
			}
		}
		theobj = theobj->next;
		iobj++;
	}
	for (inum = 0; inum < nnum; inum++) {
		if (domedian) {
			fout[inum + 1] = (double) fmedian(fmed[inum], count);
			free(fmed[inum]);
		} else {
			fout[inum + 1] /= count;
		}
	}
}
