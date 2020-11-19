/*
 * grid.c
 */

#include <stdio.h>
#include <math.h>
#include "../../catlib/cat.h"
#include "../../utils/error.h"
#include "grid.h"

static 	int	ng[2];
static	double	dx[2], x0[2];

#define	NG_MAX 100


void	setgridsize(double *xmin, double *xmax, double  d)
{
	int	j;

	for (j = 0; j < 2; j++) {
		/* first we stretch xmax a little */
       		xmax[j] += 0.01 * (xmax[j] - xmin[j]);
		ng[j] = (int) ceil((xmax[j] - xmin[j]) / d);
		dx[j] = d;
		x0[j] = xmin[j];
		if (ng[j] > NG_MAX) {
			ng[j] = NG_MAX;
			dx[j] = (xmax[j] - xmin[j]) / NG_MAX;
		}
	}
}

void	allocgrid(object ****gridptr)
{
	int ix;

	*gridptr = (object ***) calloc(ng[0], sizeof(object **));
	for (ix = 0; ix < ng[0]; ix++) {
		(*gridptr)[ix] = (object **) calloc(ng[1], sizeof(object *));
	}
}

int	getgridcoords(double *x, int *ix, int *iy)
{
	
	*ix = (int) floor((x[0] - x0[0]) / dx[0]);
	*iy = (int) floor((x[1] - x0[1]) / dx[1]);
	if (*ix < 0 || *ix >= ng[0] || *iy < 0 || *iy >= ng[1])
		return(0);
	else
		return (1);
}

void	getneighbours(double *x, object ***grid, int *nneighbours, object **neighbourobj)
{
	int	n = 0, ix0, iy0, ix, iy;
	object	*theobject;

	getgridcoords(x, &ix0, &iy0);
	for (ix = ix0 - 1; ix <= ix0 + 1; ix++) {
		if (ix < 0 || ix >= ng[0])
			continue;
		for (iy = iy0 - 1; iy <= iy0 + 1; iy++) {
			if (iy < 0 || iy >= ng[1])
				continue;
			theobject = grid[ix][iy];
			while (theobject) {
				neighbourobj[n++] = theobject;
				if (n == MAX_NEIGHBOURS)
					error_exit("mergecats: too many neighbours\n");
				theobject = theobject->next;
			}
		}
	}
	*nneighbours = n;
}



void	getobjects(object ***grid, int *nobjs, object **theobjlist)
{
	int	n = 0, ix, iy;
	object	*theobject;

	for (ix = 0; ix < ng[0]; ix++) {
		for (iy = 0; iy < ng[1]; iy++) {
			theobject = grid[ix][iy];
			while (theobject) {
				theobjlist[n++] = theobject;
				/* if (n == MAX_NEIGHBOURS)
					error_exit("mergecats: too many objects\n"); */
				theobject = theobject->next;
			}
		}
	}
	*nobjs = n;
}

