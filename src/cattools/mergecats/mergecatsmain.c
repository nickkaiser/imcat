/*
 * mergecatsmain.c
 */

#define usage "\n\n\
NAME\n\
	mergecats - merge catalogues of objects by position\n\
\n\
SYNOPSIS\n\
	mergecats [options...] d a.cat b.cat ....\n\
\n\
DESCRIPTION\n\
	'mergecats' reads N catalogues of objects and outputs a single\n\
	merged catalogue of objects whose positions match to\n\
	within tolerance d.\n\
\n\
	We first read all the catalogues\n\
	and then for each object in turn construct an N-tuplet\n\
	consisting of it and any neighbours which meet the\n\
	positional tolerance criterion.  We then rank the N-tuplets\n\
	in order of quality of match (an N-tuplet\n\
	with all slots filled ranks higher than one with one\n\
	empty slot etc., otherwise rank is the sum of the\n\
	N (N - 1) / 2 separations).  We then output the\n\
	N-tuplets (as objects with same named items as\n\
	the input catalogue but where each item is a N-vector\n\
	of the input values) in order of decreasing rank, but\n\
	only using the objects which were not contained in a previously\n\
	output N-tuplet.\n\
\n\
	By default, mergecats will only output complete ntuplets (i.e\n\
	those with detections in all input catalogues).\n\
\n\
	The idea here is that if one has three input catalogues\n\
	containing positionally coincidental objects B,V,I say, plus\n\
	an extra nearby neighbour N detected in B only,\n\
	then the algorithm will construct four triplets BVI, VIB, IBV and\n\
	NVI, it will then output whichever of the first 3 triplets\n\
	is tightest and then output an extra object N-- with\n\
	two empty slots.\n\
\n\
	For efficiency we read the objects from each catalogue into\n\
	a checkerboard grid of null terminated linked lists of objects.\n\
	Options are:\n\
\n\
	-x xname	Supply name for the 2-vector spatial coord ('x')\n\
\n\
	-n nmin		Output only objects with >= nmin detections.\n\
			With this it may be useful to use -m option:\n\
\n\
	-N nmax		Output only objects with <= nmax detections.\n\
			With this it may be useful to use -m option:\n\
\n\
	-m		Prepend the output object items with a mask which\n\
		 	is a binary representation of the detections. E.g\n\
			mask = '10010' indicates a detection in the zeroth\n\
			and third catalogues of a five catalogue merge.\n\
			Leading zeros are not printed.\n\
\n\
	-M mask		Output only objects which match the specified mask.\n\
\n\
	-s		Prepend the output cat with a column containing\n\
			the 'size' of the object (sum of the N (N - 1) / 2\n\
			separations.\n\
\n\
	-d		Prepend the output cat with a column containing\n\
			the number of detections.\n\
\n\
	-e		Exclude zero separation ntuplets\n\
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

#define BIG_POS	1.e100
#define BIG_NEG -1.e100
typedef struct ntuplet {
	int	mask;
	int	nobjs;
	double	size;
	object	**obj;
	struct	ntuplet *next;
} ntuplet;


int	identical(cathead *cat1, cathead *cat2);
void	fan(int j);
int	valid(object **obj);
double	sep2(object *obj1, object *obj2);
double	getsize(object **obj);
int	getnobjs(object **obj);
int	ntupletcmp(const void *ptr1, const void *ptr2);
void	dispose(ntuplet **ntuparray, object *opobj);

static int	N, *nneighbours, xindex, nntuplets, nobjsmin, nobjsmax,
		domask, nprepends, dosize, dondet, nozerosep, themask;
static double	d;
static object	***neighbourobj, **obj, *nullipobj;
static ntuplet	*ntupletbase;


main(int argc, char *argv[])
{
	int		arg, i, j, ix, iy, i0, dim[10], nobj;
	char		*xname, defaultxname[2] = "x", **catname;
	double		*x, xmin[2], xmax[2];
	FILE		**catf;
	cathead		**ipcat, *opcat;
	object		**ipobjbase, *ipobj, *opobj, ****grid, *nextobj;
	item		*xitem, *ipitem, *opitem;
	ntuplet		*ntup, **ntuparray;
	int		*cattype;

	/* defaults */
	xname = defaultxname;
	nobjsmin = -1;
	nobjsmax = -1;
	nprepends = 0;
	domask = dosize = dondet = 0;
	nozerosep = 0;
	themask = 0;

	/* parse args */
	if (argc < 3) {
		error_exit(usage);
	}
	arg = 1;
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			break;
		}
		switch (argv[arg++][1]) {
			case 'x':
				xname = argv[arg++];
				break;
			case 'n':
				if (1 != sscanf(argv[arg++], "%d", &nobjsmin));
				break;
			case 'N':
				if (1 != sscanf(argv[arg++], "%d", &nobjsmax));
				nobjsmin = 0;
				break;
			case 'M':
				if (1 != sscanf(argv[arg++], "%d", &themask));
 				break;
			case 'm':
				domask = 1;
				break;
			case 's':
				dosize = 1;
				break;
			case 'd':
				dondet = 1;
				break;
			case 'e':
				nozerosep = 1;
				break;
			default:
				error_exit(usage);
				break;
		}	
	}

	if (arg > argc - 3)
		error_exit(usage);

	/* get the tolerance */
	if (1 != sscanf(argv[arg++], "%lf", &d))
		error_exit(usage);

	/* N = number of cats */
	N = argc - arg;
	if (nobjsmin < 0)
		nobjsmin = N;
	if (nobjsmax < 0)
		nobjsmax = N;
	
	/* open the cats */
	catf = (FILE **) calloc(N, sizeof(FILE *));
	catname = (char **) calloc(N, sizeof(char *));
	for (i = 0; i < N; i++) {
		if (!(catf[i] = fopen(catname[i] = argv[arg++], "r"))) {
			fprintf(stderr, "mergecats: Can't open %s\n", argv[--arg]);
			exit(-1);
		}
	}


	/* read the cat headers */
	ipcat = (cathead **) calloc(N, sizeof(cathead *));
	cattype = (int *) calloc(N, sizeof(int));
	for (i = 0; i < N; i++){
		setcatipf(catf[i]);
		ipcat[i]  = readcathead();
		getcatipfiletype(&(cattype[i]));
		if (i > 0) {
			if (!identical(ipcat[i], ipcat[0])) {
				error_exit("mergecats: input cats must contain same object items\n");
			}
		}
	}

	/* check that it's a 2-vector */
	xitem = getobjectitem(xname, ipcat[0]);
	if ((xitem->itype != NUM_TYPE) || (xitem->ndim != 1) || ((xitem->dim)[0] != 2))
		error_exit("mergecats: position variable must be numerical 2-vector\n");

	/* now we create the output cathead */
	opcat = (cathead *) calloc(1, sizeof(cathead));
	copyheaderinfo(opcat, ipcat[0]);


	/* add comment, headeritems */
	addargscomment(argc, argv, opcat);
	opitem = newitem("merged_cat_names", TEXT_TYPE, 1, N);
	allocitemcontents(opitem, &(opitem->addr), 0);
	for (i = 0; i < N; i++) {
		((char **) (opitem->addr))[i] = catname[i];
	}
	installitem(opitem, &(opcat->headeritembase));
	opitem = newitem("merge_tolerance", NUM_TYPE, 1, 1);
	allocitemcontents(opitem, &(opitem->addr), 0);
	*((double *) (opitem->addr)) = d;
	installitem(opitem, &(opcat->headeritembase));

	if (domask) {
		opitem = newitem("mask", NUM_TYPE, 1, 1);
		addobjectitem(opitem, opcat);
		nprepends++;
	}
	if (dosize) {
		opitem = newitem("size", NUM_TYPE, 1, 1);
		addobjectitem(opitem, opcat);
		nprepends++;
	}
	if (dondet) {
		opitem = newitem("ndet", NUM_TYPE, 1, 1);
		addobjectitem(opitem, opcat);
		nprepends++;
	}

	/* make 2-vector obj items */
	ipitem = ipcat[0]->objectitembase;
	while (ipitem) {
		if (ipitem->ndim > 10) {
			error_exit("mergecats: I do not support > 10 dimensional arrays\n");
		}
		for (i = 0; i < ipitem->ndim; i++) {
			dim[i] = (ipitem->dim)[i];
		}
		opitem = newitem(ipitem->name, ipitem->itype, ipitem->ndim + 1, N,
			dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], dim[8], dim[9]);
		addobjectitem(opitem, opcat);
		ipitem = ipitem->next;
	}

	/* and write it out */
	setcatopfiletype(BINARY_FILE_TYPE);
	writecathead(opcat);
	/* and create the opobj */
	opobj = newobject(opcat);
	allocobjectcontents(opobj);


	/* now we read the objects into null terminated linked lists */
	/* and get xmax, xmin */
	ipobj = newobject(ipcat[0]);
	xindex = getobjectitemindex(xname, ipobj);
	xmin[0] = xmin[1] = BIG_POS;
	xmax[0] = xmax[1] = BIG_NEG;
	ipobjbase = (object **) calloc(N, sizeof(object *));
	for (i = 0; i < N; i++) {
		setcatipf(catf[i]);
		setcatipfiletype(cattype[i]);
		while (1) {
			ipobj = newobject(ipcat[0]);
			allocobjectcontents(ipobj);
			if (!readobject(ipobj))
				break;
			ipobj->next = ipobjbase[i];
			ipobjbase[i] = ipobj;
			x = (double *) ((ipobj->addrlist)[xindex]);
			for (j = 0; j < 2; j++) {
				xmin[j] = (x[j] < xmin[j] ? x[j] : xmin[j]);
				xmax[j] = (x[j] > xmax[j] ? x[j] : xmax[j]);
			}
		}
	}

	/* now we figure the size of the grid */
	setgridsize(xmin, xmax, d);

	/* create the grids of object ptrs */
	grid = (object ****) calloc(N, sizeof(object ***));
	for (i = 0; i < N; i++) {
		allocgrid(&(grid[i]));
	}

	/* install the particles in the grid */
	for(i = 0; i < N; i++) {
		ipobj = ipobjbase[i];
		while (ipobj) {
			if (!getgridcoords((double *)((ipobj->addrlist)[xindex]), &ix, &iy)) {
				ipobj = ipobj->next;
				continue;
			}
			nextobj = ipobj->next;
			ipobj->next = grid[i][ix][iy];
			grid[i][ix][iy] = ipobj;
			ipobj = nextobj;
		}
	}

	/* now make null terminated linked list of all ntuplets */
	nntuplets = 0;
	obj = (object **) calloc(N, sizeof(object *));
	nneighbours = (int *) calloc(N, sizeof(int));
	neighbourobj = (object ***) calloc(N, sizeof(object **));
	for (i = 0; i < N; i++)
		neighbourobj[i] = (object **) calloc(MAX_NEIGHBOURS, sizeof(object *));
	for (i = 0; i < N; i++) {
		getobjects(grid[i], &nobj, neighbourobj[i]);
		nneighbours[i] = -1;
		for (i0 = 0; i0 < nobj; i0++) {
			obj[i] = neighbourobj[i][i0];
			x = (double *)(((obj[i])->addrlist)[xindex]);
			for (j = 0; j < N; j++) {
				if (j == i)
					continue;
				getneighbours(x, grid[j], &(nneighbours[j]), neighbourobj[j]);
			}
			fan(0);		/* recursive function */
		}
	}

	/* transfer the ntuplets to a big array */
	ntuparray = (ntuplet **) calloc(nntuplets, sizeof(ntuplet *));
	ntup = ntupletbase;
	i = 0;
	while (ntup) {
		ntuparray[i++] = ntup;
		ntup = ntup->next;
	}

	/* and sort them */
	qsort((char *) ntuparray, nntuplets, sizeof(void *), ntupletcmp);

	/* and dispose of them */
	nullipobj = newobject(ipcat[0]);
	allocobjectcontents(nullipobj);
	dispose(ntuparray, opobj);
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



void	fan(int j)
{
	int	i, jj, mult;
	ntuplet *thentuplet;

	if (nneighbours[j] < 0) {
		if (j < (N - 1)) {
			fan(j + 1);
      		} else {
			return;
		}
	}
	for (i = -1; i < nneighbours[j]; i++) {
		if (i < 0)
			obj[j] = NULL;
		else
			obj[j] = neighbourobj[j][i];
		if (j == (N - 1) || ((j == (N - 2)) && (nneighbours[N - 1] < 0))) {
			if (valid(obj)) {
				thentuplet = (ntuplet *) calloc(1, sizeof(ntuplet));
				thentuplet->obj = (object **) calloc(N, sizeof(object *));
				mult = 1;
				for (jj = N - 1; jj >= 0; jj--) {
					thentuplet->obj[jj] = obj[jj];
					if (obj[jj])
						thentuplet->mask += mult;
					mult *= 10;
				}
				thentuplet->nobjs = getnobjs(obj);
				thentuplet->size = getsize(obj);
				thentuplet->next = ntupletbase;
				ntupletbase = thentuplet;
				nntuplets++;
			}
		} else {
			fan(j + 1);
		}
	}
}



int	valid(object **obj)
{
	int	i, j;
	double	rr;

	for (i = 0; i < N; i++) {
		if (!(obj[i]))
			continue;
		for (j = 0; j < N; j++) {
			if ((j <= i) || !(obj[j]))
				continue;
			rr = sep2(obj[i], obj[j]);
			if ((rr == 0.0) && nozerosep)
				return (0);
			if ((d * d) < rr)
				return (0);
		}
	}

	return (1);
}



double	sep2(object *obj1, object *obj2)
{
	double	*x1, *x2;
	double	dx, dy;

	x1 = (double *) ((obj1->addrlist)[xindex]);
	x2 = (double *) ((obj2->addrlist)[xindex]);
	dx = x2[0] - x1[0];
	dy = x2[1] - x1[1];
	return (dx * dx + dy * dy);
}


double	getsize(object **obj)
{
	int	i, j;
	double	size = 0.0;

	for (i = 0; i < N; i++) {
		if (!(obj[i]))
			continue;
		for (j = i + 1; j < N; j++) {
			if (!(obj[j]))
				continue;
			size += sqrt(sep2(obj[i], obj[j]));
		}
	}
	return (size);
}


int	getnobjs(object **obj)
{
	int 	i, nobjs = 0;

	for (i = 0; i < N; i++) {
		if (obj[i]) {
			nobjs++;
		}
	}
	return (nobjs);
}


#define query(a,b) ((a) > (b) ? 1 : ((a) < (b) ? -1 : 0))
int	ntupletcmp(const void *ptr1, const void *ptr2)
{
	ntuplet **ntup1, **ntup2;
	int	nobjstest;

	ntup1 = (ntuplet **) ptr1;
	ntup2 = (ntuplet **) ptr2;
	nobjstest = query((*ntup2)->nobjs, (*ntup1)->nobjs);
	return ((nobjstest ? nobjstest : query((*ntup1)->size, (*ntup2)->size))); 
}
#undef query


void	dispose(ntuplet **ntuparray, object *opobj)
{
	object  *ipobj, *theobj;
	int	i0, i, j, badun;
	item	*theitem;

	for (i0 = 0; i0 < nntuplets; i0++) {
		if ((ntuparray[i0])->nobjs < nobjsmin)
			exit(0);
		badun = 0;
		for (j = 0; j < N; j++) {
			ipobj = ((ntuparray[i0])->obj)[j];
			if (ipobj) {
				if (!(ipobj->addrlist)) {
					badun = 1;
				}
			}
		}
		if (badun) {
			continue;
		}
		i = 0;
		if (domask) {
			*((double *) (opobj->addrlist)[i]) = (ntuparray[i0])->mask;
			i++;
		}
		if (dosize) {
			*((double *) (opobj->addrlist)[i]) = (ntuparray[i0])->size;
			i++;
		}
		if (dondet) {
			*((double *) (opobj->addrlist)[i]) = (ntuparray[i0])->nobjs;
			i++;
		}
		for (i = 0; i < (opobj->nitems - nprepends); i++) {
			theitem = ((opobj->cathead)->itemlist)[i + nprepends];
			for (j = 0; j < N; j++) {
				ipobj = ((ntuparray[i0])->obj)[j];
				theobj = (ipobj ? ipobj : nullipobj);
				((void ***)(opobj->addrlist))[i + nprepends][j] = (theobj->addrlist)[i];
			}
		}
		if ((ntuparray[i0])->nobjs <= nobjsmax) {
			if (!themask || ((ntuparray[i0])->mask == themask)) {
				writeobject(opobj);
			}
		}
		for (j = 0; j < N; j++) {
			ipobj = ((ntuparray[i0])->obj)[j];
			if (ipobj)
				ipobj->addrlist = NULL;
		}
		
	}
	exit(0);
}
