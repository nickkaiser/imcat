 /*
 * pair.c
 */

#define usage "\n\n\
NAME\n\
	getisolatedobjects - find isolated objects\n\
\n\
SYNOPSIS\n\
	getisolatedobjects [options....] dmax referencecat\n\
		-u		# print this man page\n\
		-x xname 	# specify name for coords (x)\n\
		-e		# don't exclude zero separation objects\n\
\n\
DESCRIPTION\n\
        'getisolatedobjects' reads a source catalogue from stdin, a\n\
	reference catalogue from 'referencecat' and writes to stdout\n\
	a catalogue containing objects from the source cat which\n\
	have no reference neighbour with |separation| < d.\n\
	It works by installing the objects from the reference cat into\n\
	a grid of linked lists, and it then processes the objects\n\
	from the source cat sequentially.\n\
	The -e option is useful if one wants to find the objects\n\
	in a cat which have no close neighbours in the same catalogue.\n\
	To do this simply use the source cat as reference cat.\n\
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

double	sep2(object *ipobj, int ipxindex, object *refobj, int refxindex);

#define BIG_POS	1.e100
#define BIG_NEG -1.e100

main(int argc, char *argv[])
{
	int		arg = 1, dim[10], i, j, index, ipxindex, refxindex, ix, iy, nneighbours;
	int		refcatfiletype, ipcatfiletype;
	cathead		*refcat, *thecat;
	FILE		*refcatf, *ipcatf = stdin;
	object		*ipobj, *refobj, *refobjbase = NULL, *obj1, *obj2, ***grid, *nextobj;
	object		*neighbourobject[MAX_NEIGHBOURS];
	item		*ipxitem, *refxitem;
	double		dmax, *x, xmin[2], xmax[2], sepsquared;
	char		*xname, defaultxname[2] = "x";
	int		excludezerosep, isolated;

	/* defaults */
	xname = defaultxname;
	excludezerosep = 1;

	/* parse the args */
	while (arg < argc) {
		if (argv[arg][0] != '-')
			break;
		switch (argv[arg++][1]) {
			case 'x':
				xname = argv[arg++];
				break;
			case 'e':
				excludezerosep = 0;
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	/* read the max separation */
	if (1 != sscanf(argv[arg++], "%lf", &dmax))
		error_exit(usage);

	/* try to open the cats */
	if (!(refcatf = fopen(argv[arg], "r"))) {
		fprintf(stderr, "getisolatedobjects: cannot open %s\n", argv[arg]);
		exit(-1);
	}

	/* read the cat headers */
	setcatipf(refcatf);
	refcat = readcathead();
	getcatipfiletype(&refcatfiletype);
	setcatipf(ipcatf);
	thecat = readcathead();
	getcatipfiletype(&ipcatfiletype);

	/* check that they both contain 2-vector coord */
	ipxitem = getobjectitem(xname, thecat);
	if ((ipxitem->itype != NUM_TYPE) || (ipxitem->ndim != 1) || ((ipxitem->dim)[0] != 2))
		error_exit("getisolatedobject: source cat position variable must be numerical 2-vector\n");
	refxitem = getobjectitem(xname, refcat);
	if ((refxitem->itype != NUM_TYPE) || (refxitem->ndim != 1) || ((refxitem->dim)[0] != 2))
		error_exit("getisolatedobject: reference cat position variable must be numerical 2-vector\n");

	/* create the output cathead
	opcat = (cathead *) calloc(1, sizeof(cathead));
	copyheaderinfo(opcat, ipcat);
 */
	/* add comment */
	addargscomment(argc, argv, thecat);

	/* and write output cathead */
	setcatopfiletype(BINARY_FILE_TYPE);
	writecathead(thecat);


	/* read the refcat into a linked list */
	setcatipf(refcatf);
	setcatipfiletype(refcatfiletype);
	refobj = newobject(refcat);
	refxindex = getobjectitemindex(xname, refobj);
	xmin[0] = xmin[1] = BIG_POS;
	xmax[0] = xmax[1] = BIG_NEG;
	while (1) {
		refobj = newobject(refcat);
		allocobjectcontents(refobj);
		if (!readobject(refobj))
			break;
		refobj->next = refobjbase;
		refobjbase = refobj;
		x = (double *) ((refobj->addrlist)[refxindex]);
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
	refobj = refobjbase;
	while (refobj) {
		if (!getgridcoords((double *)((refobj->addrlist)[refxindex]), &ix, &iy)) {
			refobj = refobj->next;
			continue;
		}
		nextobj = refobj->next;
		refobj->next = grid[ix][iy];
		grid[ix][iy] = refobj;
		refobj = nextobj;
	}

	/* loop over objects from source cat */
	setcatipf(ipcatf);
	setcatipfiletype(ipcatfiletype);
	ipobj = newobject(thecat);
	ipxindex = getobjectitemindex(xname, ipobj);
	allocobjectcontents(ipobj);
	while (readobject(ipobj)) {
		x = (double *) ((ipobj->addrlist)[ipxindex]);
		getneighbours(x, grid, &nneighbours, neighbourobject);
		isolated = 1;
		for (i = 0; i < nneighbours; i++) {
			refobj = neighbourobject[i];
			sepsquared = sep2(ipobj, ipxindex, refobj, refxindex);
			if (sepsquared < (dmax * dmax)) {
				if (sepsquared > 0.0 || excludezerosep) {
					isolated = 0;
					break;
				}
			}
		}
		if (isolated)
			writeobject(ipobj);
	}
	exit(0);
}







double	sep2(object *ipobj, int ipxindex, object *refobj, int refxindex)
{
	double	*refx, *ipx;
	double	dx, dy;

	refx = (double *) ((refobj->addrlist)[refxindex]);
	ipx = (double *) ((ipobj->addrlist)[ipxindex]);
	dx = ipx[0] - refx[0];
	dy = ipx[1] - refx[1];
	return (dx * dx + dy * dy);
}

