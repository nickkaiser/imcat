/*
 * cleancat.c
 */

#define usage "\n\n\
NAME\n\
	cleancat - remove objects with brighter close neighbours\n\
\n\
SYNOPSIS\n\
	cleancat [options....] dmax\n\
		-u		# print this man page\n\
		-x xname 	# specify name for coords (x)\n\
		-m magname 	# specify name for magnitude (mag)\n\
\n\
DESCRIPTION\n\
        'cleancat' reads a catalogue from stdin and writes to stdout\n\
	a catalogue containing only those objects from the input\n\
	catalogue which have no brighter neighbours within distance\n\
	dmax.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n"


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "../../utils/error.h"
#include "../../catlib/cat.h"
#include "grid.h"

double	sep2(object *obj1, object *obj2, int xindex);

#define BIG_POS	1.e100
#define BIG_NEG -1.e100

main(int argc, char *argv[])
{
	int		arg = 1, dim[10], i, j, index, xindex, magindex, ix, iy, nneighbours, nobjects, nn;
	int		ipfiletype, objectisclean;
	cathead		*cathead;
	object		*ipobj, *ipobjbase = NULL, *opobj, *obj1, *obj2, ***grid, *nextobj, **objlist;
	object		*neighbourobject[MAX_NEIGHBOURS];
	item		*ipitem, *opitem, *xitem, *magitem;
	void		**opaddr;
	double		d, dmax, *x, xmin[2], xmax[2], *mag1, *mag2;
	char		*xname, defaultxname[2] = "x", *magname, defaultmagname[4] = "mag";

	/* defaults */
	xname = defaultxname;
	magname = defaultmagname;

	/* parse the args */
	while (arg < argc) {
		if (argv[arg][0] != '-')
			break;
		switch (argv[arg++][1]) {
			case 'x':
				xname = argv[arg++];
				break;
			case 'm':
				magname = argv[arg++];
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	/* read the max separation */
	if (1 != sscanf(argv[arg++], "%lf", &dmax))
		error_exit(usage);


	/* read the cat header */
	cathead = readcathead();
 	getcatipfiletype(&ipfiletype);
	
	/* check that it contains a 2-vector coord and magnitude*/
	xitem = getobjectitem(xname, cathead);
	if ((xitem->itype != NUM_TYPE) || (xitem->ndim != 1) || ((xitem->dim)[0] != 2))
		error_exit("cleancat: position variable must be numerical 2-vector\n");
	magitem = getobjectitem(magname, cathead);
	if ((magitem->itype != NUM_TYPE) || (magitem->ndim != 1) || ((magitem->dim)[0] != 1))
		error_exit("cleancat: magnitude variable must be a scalar\n");

	/* add comment */
	addargscomment(argc, argv, cathead);

	/* and write output cathead */
	setcatopfiletype(ipfiletype);
	writecathead(cathead);

	/* read the cat into a linked list */
	nobjects = 0;
	ipobj = newobject(cathead);
	xindex = getobjectitemindex(xname, ipobj);
	magindex = getobjectitemindex(magname, ipobj);
	xmin[0] = xmin[1] = BIG_POS;
	xmax[0] = xmax[1] = BIG_NEG;
	while (1) {
		ipobj = newobject(cathead);
		allocobjectcontents(ipobj);
		if (!readobject(ipobj))
			break;
		nobjects++;
		ipobj->next = ipobjbase;
		ipobjbase = ipobj;
		x = (double *) ((ipobj->addrlist)[xindex]);
		for (j = 0; j < 2; j++) {
			xmin[j] = (x[j] < xmin[j] ? x[j] : xmin[j]);
			xmax[j] = (x[j] > xmax[j] ? x[j] : xmax[j]);
		}
	}

	/* now we figure the size of the grid */
	setgridsize(xmin, xmax, dmax);

	/* create the grid of object ptrs */
	allocgrid(&grid);

	/* install the particles in the grid */
	ipobj = ipobjbase;
	while (ipobj) {
		if (!getgridcoords((double *)((ipobj->addrlist)[xindex]), &ix, &iy)) {
			ipobj = ipobj->next;
			continue;
		}
		nextobj = ipobj->next;
		ipobj->next = grid[ix][iy];
		grid[ix][iy] = ipobj;
		ipobj = nextobj;
	}

	/* get array of pointers to the objects */
	objlist = (object **) calloc(nobjects, sizeof(double));
	getobjects(grid, &nn, objlist);

	/* loop over objects in list and get neighbours*/
	for (j = 0; j < nobjects; j++) {
		obj1 = objlist[j];
		x = (double *) ((obj1->addrlist)[xindex]);
		getneighbours(x, grid, &nneighbours, neighbourobject);
		objectisclean = 1;
		for (i = 0; i < nneighbours; i++) {
			obj2 = neighbourobject[i];
			d = sqrt(sep2(obj1, obj2, xindex));
			if (d < dmax && d > 0.0) {
				mag1 = (double *) ((obj1->addrlist)[magindex]);
				mag2 = (double *) ((obj2->addrlist)[magindex]);
				if (*mag2 < *mag1) {
					objectisclean = 0;
					break;
				}
			}
		}
		if (objectisclean) {
			writeobject(obj1);
		}
	}
	exit(0);
}








double	sep2(object *obj1, object *obj2, int xindex)
{
	double	*x1, *x2;
	double	dx, dy;

	x1 = (double *) ((obj1->addrlist)[xindex]);
	x2 = (double *) ((obj2->addrlist)[xindex]);
	dx = x2[0] - x1[0];
	dy = x2[1] - x1[1];
	return (dx * dx + dy * dy);
}

