 /*
 * findtrailobjects.c
 */

#define usage "\n\n\
NAME\n\
	findtrailobjects - find objects in very narrow satellite trails etc\n\
\n\
SYNOPSIS\n\
	findtrailobjects [options....] dmax\n\
		-u		# print this man page\n\
		-x xname 	# specify name for coords (x)\n\
		-n nphi		# number of bins for histogram (512)\n\
\n\
DESCRIPTION\n\
        'findtrailobjects' reads a catalogue from stdin which must contain\n\
	a 2 vector 'x' (substitute some other name with -x option).\n\
	For each object it computes a histogram of angles to neighbours\n\
	within distance dmax, where angle lies in domain 0 - PI.\n\
	It outputs a lc binary format catalogue containing\n\
		x		# the position\n\
		nmax		# the count in the highest histogram bin\n\
		phi		# the angle of the highest histogram bin\n\
	With judicious choice of parameters and filtering, this will\n\
	detect satellite trails, diffraction spikes around stars etc.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n"


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "../../utils/error.h"
#include "../../utils/args.h"
#include "../../catlib/cat.h"
#include "../../imlib/fits.h"
#include "grid.h"

#define	NOPBUFF	4

main(int argc, char *argv[])
{
	char	*flag, *xname, defaultxname[2] = "x", lccommand[128], argstring[512];
	double	*x, *x1, *x2, xmin[2], xmax[2], dmax, d, dphi, dx, dy, phi, phi0, phimax;
	cathead	*ipcat;
	object	*theobject, *nextobject, *baseobject, ***grid;
	object	**objlist, *neighbourobject[MAX_NEIGHBOURS], *obj1, *obj2;
	item	*xitem;
	int	xindex, i, j, ix, iy, nobjects, nn, nneighbours;
	int	*n, iphi, nphi, nmax;
	double	opbuff[NOPBUFF];

	/* defaults */
	xname = defaultxname;
	nphi = 512;

	/* parse the args */
	argsinit(argc, argv, usage);
	while (nextargtype() == FLAG_ARG) {
		flag = getflag();
		switch (flag[0]) {
			case 'x':
				xname = getargs();
				break;
			case 'n':
				nphi = getargi();
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}
	dmax = getargd();

	/* set up the histogram */
	dphi = M_PI / nphi;
	n = (int *) calloc(nphi, sizeof(int));

	/* read the cat header */
	ipcat = readcathead();

	/* check that it contains the 2-vector coord */
	xitem = getobjectitem(xname, ipcat);
	if ((xitem->itype != NUM_TYPE) || (xitem->ndim != 1) || ((xitem->dim)[0] != 2)) {
		error_exit("findtrailobjects: position variable must be numerical 2-vector\n");
	}

	/* ouput the catalogue header */
	argsToString(argc, argv, argstring);
	sprintf(lccommand, "lc -C -b -N '1 2 %s' -n nmax -n phi -x -a '%s' < /dev/null", xname, argstring);
	system(lccommand);

	/* read the cat into a linked list */
        nobjects = 0;	
	theobject = newobject(ipcat);
	xindex = getobjectitemindex(xname, theobject);
	xmin[0] = xmin[1] = FLT_MAX;
	xmax[0] = xmax[1] = -FLT_MAX;
	baseobject = NULL;
	while (1) {
		theobject = newobject(ipcat);
		allocobjectcontents(theobject);
		if (!readobject(theobject))
			break;
		nobjects++;
		theobject->next = baseobject;
		baseobject = theobject;
		x = (double *) ((theobject->addrlist)[xindex]);
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
	theobject = baseobject;
	while (theobject) {
		if (!getgridcoords((double *)((theobject->addrlist)[xindex]), &ix, &iy)) {
			theobject = theobject->next;
			continue;
		}
		nextobject = theobject->next;
		theobject->next = grid[ix][iy];
		grid[ix][iy] = theobject;
		theobject = nextobject;
	}

        /* get array of pointers to all the objects */
        objlist = (object **) calloc(nobjects, sizeof(double));
        getobjects(grid, &nn, objlist);

        /* loop over pairs with separation less than d and compute histograms */
	for (j = 0; j < nobjects; j++) {
 		/* reset the histogram */
		for (iphi = 0; iphi < nphi; iphi++) {
			n[iphi] = 0;
		}
		obj1 = objlist[j];
		x1 = (double *) ((obj1->addrlist)[xindex]);
		getneighbours(x1, grid, &nneighbours, neighbourobject);
		for (i = 0; i < nneighbours; i++) {
			obj2 = neighbourobject[i];
			x2 = (double *) ((obj2->addrlist)[xindex]);
        		dx = x2[0] - x1[0];
        		dy = x2[1] - x1[1];
			if (dmax > sqrt(dx * dx + dy * dy)) {
				/* get the angle on -PI < phi < PI domain */
				phi0 = phi = atan2(dy, dx);
				/* get the angle on 0 < phi0 < PI domain */
				if (phi0 < 0.0) {
					phi0 += M_PI;
				}
				/* shift by 0.5 * dphi so 1st bin straddles phi0 = 0 */
				phi0 += 0.5 * dphi;
				if (phi0 > M_PI) {
					phi0 -= M_PI;
				}
				iphi = floor(phi0 /= dphi);
				if (iphi >= 0 && iphi < nphi) {
					n[iphi]++;
				}
			}
		}
		/* find highest bin */
		nmax = 0;
		for (iphi = 0; iphi < nphi; iphi++) {
			if (n[iphi] > nmax) {
				nmax = n[iphi];
				phimax = dphi * (iphi + 0.5);
			}
		}
		opbuff[0] = x1[0];
		opbuff[1] = x1[1];
		opbuff[2] = (double) nmax;
		opbuff[3] = phimax;
		fwrite(opbuff, sizeof(double), NOPBUFF, stdout);
	}

	exit(0);
}



