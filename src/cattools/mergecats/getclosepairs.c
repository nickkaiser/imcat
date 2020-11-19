 /*
 * pair.c
 */

#define usage "\n\n\
NAME\n\
	getclosepairs - make catalogue of close pairs from 2 input cats\n\
\n\
SYNOPSIS\n\
	getclosepairs [options....] dmax cat1 cat2\n\
		-u		# print this man page\n\
		-x xname 	# specify name for coords (x)\n\
\n\
DESCRIPTION\n\
        'getclosepairs' reads two catalogues cat1, cat2 and writes to stdout\n\
	a catalogue containing pairs with |separation| < dmax \n\
	Output cat has object items with the same names as in the input cats\n\
	(which must have identical object items)\n\
	but each object is a 2-vector formed from a pair of\n\
	input objects.\n\
	It works by installing the objects from the first cat into\n\
	a grid of linked lists, and it then reads the objects\n\
	from the second cat and outputs a pair-object for\n\
	each pair with separation < dmax.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n"


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "../../utils/error.h"
#include "../../catlib/cat.h"
#include "grid.h"

int	identical(cathead *cat1, cathead *cat2);
double	sep2(object *obj1, object *obj2, int xindex);

#define BIG_POS	1.e100
#define BIG_NEG -1.e100

main(int argc, char *argv[])
{
	int		arg = 1, dim[10], i, j, index, xindex, ix, iy, nneighbours;
	int		ipfiletype1, ipfiletype2;
	cathead		*cat1, *cat2, *opcat;
	FILE		*catf1, *catf2;
	object		*ipobj1, *ipobj2, *ipobjbase1 = NULL, *opobj, *obj1, *obj2, ***grid, *nextobj;
	object		*neighbourobject[MAX_NEIGHBOURS];
	item		*ipitem1, *ipitem2, *opitem, *xitem;
	void		**opaddr;
	double		dmax, *x, xmin[2], xmax[2];
	char		*xname, defaultxname[2] = "x";

	/* defaults */
	xname = defaultxname;

	/* parse the args */
	while (arg < argc) {
		if (argv[arg][0] != '-')
			break;
		switch (argv[arg++][1]) {
			case 'x':
				xname = argv[arg++];
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	/* read the max spearation */
	if (1 != sscanf(argv[arg++], "%lf", &dmax))
		error_exit(usage);

	/* try to open the cats */
	if (!(catf1 = fopen(argv[arg], "r"))) {
		fprintf(stderr, "getclosepairs: cannot open %s\n", argv[arg]);
		exit(-1);
	}
	arg++;
	if (!(catf2 = fopen(argv[arg], "r"))) {
		fprintf(stderr, "getclosepairs: cannot open %s\n", argv[arg]);
		exit(-1);
	}

	/* read the cat headers */
	setcatipf(catf1);
	cat1 = readcathead();
 	getcatipfiletype(&ipfiletype1);
	setcatipf(catf2);
	cat2 = readcathead();
	getcatipfiletype(&ipfiletype2);

	/* check that they contain identical items */
	if (!identical(cat1, cat2))
		error_exit("getclosepairs: input cats must contain the same items!\n");
	
	/* check that they contain 2-vector coord */
	xitem = getobjectitem(xname, cat1);
	if ((xitem->itype != NUM_TYPE) || (xitem->ndim != 1) || ((xitem->dim)[0] != 2))
		error_exit("mergecats: position variable must be numerical 2-vector\n");

	/* create the output cathead */
	opcat = (cathead *) calloc(1, sizeof(cathead));
	copyheaderinfo(opcat, cat1);

	/* add comment */
	addargscomment(argc, argv, opcat);


	/* create the object items */
	ipitem1 = cat1->objectitembase;
	while (ipitem1) {
		allocitemcontents(ipitem1, &(ipitem1->addr), 0);
		for (i = 0; i < ipitem1->ndim; i++) {
			dim[i] = (ipitem1->dim)[i];
		}
		opitem = newitem(ipitem1->name, ipitem1->itype, ipitem1->ndim + 1, 2,
                        dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], dim[8], dim[9]);
                addobjectitem(opitem, opcat);
               	ipitem1 = ipitem1->next;
	}

	/* and write output cathead */
	setcatopfiletype(BINARY_FILE_TYPE);
	writecathead(opcat);

	/* read the cat1 into a linked list */
	setcatipf(catf1);
	setcatipfiletype(ipfiletype1);
	ipobj1 = newobject(cat1);
	xindex = getobjectitemindex(xname, ipobj1);
	xmin[0] = xmin[1] = BIG_POS;
	xmax[0] = xmax[1] = BIG_NEG;
	while (1) {
		ipobj1 = newobject(cat1);
		allocobjectcontents(ipobj1);
		if (!readobject(ipobj1))
			break;
		ipobj1->next = ipobjbase1;
		ipobjbase1 = ipobj1;
		x = (double *) ((ipobj1->addrlist)[xindex]);
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
	ipobj1 = ipobjbase1;
	while (ipobj1) {
		if (!getgridcoords((double *)((ipobj1->addrlist)[xindex]), &ix, &iy)) {
			ipobj1 = ipobj1->next;
			continue;
		}
		nextobj = ipobj1->next;
		ipobj1->next = grid[ix][iy];
		grid[ix][iy] = ipobj1;
		ipobj1 = nextobj;
	}
	
	/* create the ouput object */
	opobj = newobject(opcat);
	allocobjectcontents(opobj);

	/* loop over objects from second cat */
	setcatipf(catf2);
	setcatipfiletype(ipfiletype2);
	ipobj2 = newobject(cat2);
	allocobjectcontents(ipobj2);
	while (readobject(ipobj2)) {
		x = (double *) ((ipobj2->addrlist)[xindex]);
		getneighbours(x, grid, &nneighbours, neighbourobject);
		for (i = 0; i < nneighbours; i++) {
			ipobj1 = neighbourobject[i];
			if (sep2(ipobj1, ipobj2, xindex) < (dmax * dmax)) {
				for (index = 0; index < opobj->nitems; index++) {
					opaddr = (void **) (opobj->addrlist)[index];
					opaddr[0] = (ipobj1->addrlist)[index];
					opaddr[1] = (ipobj2->addrlist)[index];
				}
				writeobject(opobj);
			}
		}
	}
	exit(0);
}





int	identical(cathead *cat1, cathead *cat2)
{
	item 	*objitem1, *objitem2;
	int	i, j;

	if (cat1->nobjectitems != cat2->nobjectitems) {
		return (0);
	}
	objitem1 = cat1->objectitembase;
	objitem2 = cat2->objectitembase;
	for (i = 0; i < cat1->nobjectitems; i++) {
		if (strcmp(objitem1->name, objitem2->name))
			return(0);
		if ((objitem1->itype != objitem2->itype) || (objitem1->ndim != objitem2->ndim))
			return (0);
		for (j = 0; j < objitem1->ndim; j++) {
			if ((objitem1->dim)[j] != (objitem2->dim)[j])
				return (0);
		}
		objitem1 = objitem1->next;
		objitem2 = objitem2->next;
	}

	return (1);
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

