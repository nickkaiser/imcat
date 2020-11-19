/*
 * map.c
 */

/*
 *	routines to facilitate mapping of a 'source image' onto an image
 *	on the deflector plane with the mapping supplied by a user defined
 *	function. Used to stretch images for registration or to simulate
 *	lensing.
 *
 */

#include	<stdio.h> 
#include	<stdlib.h>
#include	<math.h>
#include	<limits.h>
#include	"../utils/error.h"
#include	"fits.h"
#include	"map_private.h"
#include	"map.h"

#define	TINY	1.e-8
#define	NOCUT	0
#define	XCUT	1
#define	YCUT	2

#define       MAX(x,y) (((x) > (y)) ? (x) : (y)) 
#define       MIN(x,y) (((x) < (y)) ? (x) : (y)) 

static	int	globalmapmode = FORWARDMAPMODE;


vert	*makevertex(float x, float y)
{
	vert	*thevert;
	
/*	thevert = (vert *) calloc(1, sizeof(vert));
	if (!thevert)
		error_exit("makevertex: memory allocation failure\n");*/
	thevert = allocvert();
	thevert->x = x;
	thevert->y = y;
	return (thevert);
}



void	makering(vert **basevertex, vert **point, int npoints)
{
	int	i;
	
	for (i = 0; i < npoints - 1; i++) {
		point[i]->next = point[i+1];
		point[i+1]->prev = point[i];
	}
	point[0]->prev = point[npoints - 1];
	point[npoints - 1]->next = point[0];
	*basevertex = point[0];
}


void	decompose(vert *basevert)
/* this is kinda complicated:
 *	we want to break a general polygon up into polygons where all the
 *	(ix,iy) values are identical --- i.e they belong to the same source pixel
 *
 *	we do it recursively:
 *		we look at the polygon
 *		if all vertices belong to same source pixel
 *			we call dopolygon(basevert) and return
 *		else
 *			we break the polygon into two subpolygons by adding four new vertices:
 *				R1, L1 to close the polygon containing basevert
 *				R2, L2 to close the other subpolygon
 *			call decompose on each subpolygon
 *			
 */
{
	vert	*thevert, *r1, *r2, *l1, *l2, *R1, *R2, *L1, *L2, *next, *prev, *xmaxvert, *ymaxvert;
	int	ix0, iy0, xcrosser, upcrosser, icut;
	float	xmin, xmax, ymin, xcut, ymax, ycut;
	int	cut;
	
	/* get the min max x,y values and corresponding verts */
	thevert = xmaxvert = ymaxvert = basevert;
	xmin = xmax = basevert->x;
	ymin = ymax = basevert->y;
	while (1) {
		next = thevert->next;
		if (next->x > xmax) {
			xmax = next->x;
			xmaxvert = next;
		}
		if (next->x < xmin)
			xmin = next->x;
		if (next->y > ymax) {
			ymax = next->y;
			ymaxvert = next;
		}
		if (next->y < ymin)
			ymin = next->y;
		thevert = next;
		if (thevert->next == basevert)
			break;
	}
	
	cut = NOCUT;
	if (ceil(ymax) - floor(ymin) > 1.0) {		/* need to cut in y */
		ycut = ceil(ymax) - 1.0;
		cut = YCUT;
		thevert = ymaxvert;
		while (1) {
			next = thevert->next;
			if (next->y <= ycut) {
				r1 = thevert;
				r2 = next;
				break;
			}
			thevert = next;
		}
		thevert = ymaxvert;
		while (1) {
			prev = thevert->prev;
			if (prev->y <= ycut) {
				l1 = thevert;
				l2 = prev;
				break;
			}
			thevert = prev;
		}
	}
	if ((cut == NOCUT) && (ceil(xmax) - floor(xmin) > 1.0)) {	/* need to cut in x */
		xcut = ceil(xmax) - 1.0;
		cut = XCUT;
		thevert = xmaxvert;
		while (1) {
			next = thevert->next;
			if (next->x <= xcut) {
				r1 = thevert;
				r2 = next;
				break;
			}
			thevert = next;
		}
		thevert = xmaxvert;
		while (1) {
			prev = thevert->prev;
			if (prev->x <= xcut) {
				l1 = thevert;
				l2 = prev;
				break;
			}
			thevert = prev;
		}
	}
	
	switch (cut) {
		case NOCUT:
			dopolygon(basevert);
			return;
			break;
		case YCUT:
			makebreakpoints(r1, r2, &R1, &R2, YCUT, ycut);
			makebreakpoints(l1, l2, &L1, &L2, YCUT, ycut);
			link(r1, r2, l1, l2, R1, R2, L1, L2);
			decompose(ymaxvert);
			decompose(R2);
			return;
			break;
		case XCUT:
			makebreakpoints(r1, r2, &R1, &R2, XCUT, xcut);
			makebreakpoints(l1, l2, &L1, &L2, XCUT, xcut);
			link(r1, r2, l1, l2, R1, R2, L1, L2);
			decompose(xmaxvert);
			decompose(R2);
			return;
			break;
		default:
			error_exit("decompose: bad cut\n");
			break;
	}
}



void	dopolygon(vert *basevert)
{
	/* we now have a polygon contained within one source pixel */
	/* if it's a triangle we call dotriangle() otherwise we recurse */
	vert	*thevert, *vert1, *vert2, *prev, *next;
	int	nverts = 0;

	/* count the vertices */
	thevert = basevert;
	while (1) {
		nverts++;
		if (thevert->next == basevert)
			break;
		thevert = thevert->next;
	}
	
	if (nverts < 3)
		error_exit("dopolygon: dunno what to do with two-vertex polygon\n");
	if (nverts == 3) {
		dotriangle(basevert);
	} else {
/*		vert1 = (vert *) calloc(1, sizeof(vert));
		vert2 = (vert *) calloc(1, sizeof(vert));*/
		vert1 = allocvert();
		vert2 = allocvert();
		next = basevert->next;
		prev = basevert->prev;
		*vert1 = *next;
		*vert2 = *prev;
		(next->next)->prev = vert1;
		(prev->prev)->next = vert2;
		vert1->prev = vert2;
		vert2->next = vert1;
		next->next = prev;
		prev->prev = next;
		dotriangle(basevert);
		dopolygon(vert1);
	}
	return;
}

void	dotriangle(vert *basevert)
{
/* we should now have a triangle completely contained within a source pixel */
	float	area;
	int		i, j;
	
	area = trianglearea(basevert);
	getsourceij(basevert, &i, &j);
	addarea(i, j, area);
/*	freepolygon(basevert);*/
}



void	printpolygon(vert *basevert)
{
	vert	*thevert;
	int	nverts = 0;
	
	thevert = basevert;
	while (1) {
		nverts++;
		if (thevert->next == basevert)
			break;
		thevert = thevert->next;
	}
	fprintf(stderr, "# %2d vertex polygon\n", nverts);
	thevert = basevert;
	while (1) {
		printvertex(thevert);
		if (thevert->next == basevert)
			break;
		thevert = thevert->next;
	}
}



void	printvertex(vert * thevert)
{
	fprintf(stderr, "%10.3f %10.3f\n", 
		thevert->x, thevert->y);
}


void	smprintpolygon(vert *basevert)
{
	vert	*thevert;
	
	thevert = basevert;
	
	fprintf(stdout, "\trelocate %10.3f %10.3f\n", basevert->x, basevert->y);
	while (1) {
		thevert = thevert->next;
	fprintf(stdout, "\tdraw %10.3f %10.3f\n", thevert->x, thevert->y);
		if (thevert == basevert)
			return;
	}
}



void	freepolygon(vert *basevert)
{
	vert	*thevert, *next;
	
	thevert = basevert;
	
	while (1) {
		next = thevert->next;
		free(thevert);
		if (next == basevert);
			return;
		thevert = next;
	}
}


void	switchxy(vert *v)
{
	float	temp;
	
	temp = v->x;
	v->x = v->y;
	v->y = temp;
}



void	makebreakpoints(vert *v1, vert *v2, vert **V1, vert **V2, int dir, float cut)
{
	float	frac;
	
/*	*V1 = (vert *) calloc(1, sizeof(vert));
	*V2 = (vert *) calloc(1, sizeof(vert));*/

	*V1 = allocvert();
	*V2 = allocvert();
	
	if (dir == YCUT) {
		switchxy(v1);
		switchxy(v2);
	}
	
	(*V1)->x = (*V2)->x = cut;
	frac = (cut - v1->x) / (v2->x - v1->x);
	(*V1)->y = (*V2)->y = v1->y + frac * (v2->y - v1->y);

	if (dir == YCUT) {
		switchxy(v1);
		switchxy(v2);
		switchxy(*V1);
		switchxy(*V2);
	}
}


void	link(vert *r1, vert *r2, vert *l1, vert *l2, vert *R1, vert *R2, vert *L1, vert *L2)
{
	r1->next = R1;
	R1->prev = r1;
	R1->next = L1;
	L1->prev = R1;
	L1->next = l1;
	l1->prev = L1;
	l2->next = L2;
	L2->prev = l2;
	L2->next = R2;
	R2->prev = L2;
	R2->next = r2;
	r2->prev = R2;
}


float	trianglearea(vert *basevert)
{
	float	x1, y1, x2, y2;
	
	x1 = (basevert->next)->x - basevert->x;
	x2 = (basevert->prev)->x - basevert->x;
	y1 = (basevert->next)->y - basevert->y;
	y2 = (basevert->prev)->y - basevert->y;
	return(0.5 * fabs(x1 * y2 - x2 * y1));
}



void	getsourceij(vert *basevert, int *ix, int *iy)
{
	float	x0, x1, x2, y0, y1, y2;
	
	x0 = basevert->x;
	y0 = basevert->y;
	x1 = (basevert->next)->x;
	y1 = (basevert->next)->y;
	x2 = (basevert->prev)->x;
	y2 = (basevert->prev)->y;
	*ix = (int) floor(MIN(x0, MIN(x1, x2)));
	*iy = (int) floor(MIN(y0, MIN(y1, y2)));
}



/* globals used by map(), trianglecount()  and addarea() */
static	int 	gM1, gM2, gtargeti, gtargetj;
static	float	gareasum, gfsum, **gfsource, **gftarget;
static	long	gntriangles = 0;
static	float	garea;


/*
 * map() maps a M2 x M1 source image fsource[i][j] onto an N2 x N1 target 
 * image ftarget[i][j] with the mapping 
 *		ftarget(r) += fsource(r + d)
 * where the deflection d = (di, dj) is supplied by function deflection()
 *
 * we use the slightly confusing notation (i,j) -> (y, x)
 */
void	map(float **ftarget, int N1, int N2, float **fsource, int M1, int M2,
			int (*deflection)(float ri, float rj, float *di, float *dj))
{
	int	i, j, index, goodtriangle;
	float	di, dj;
	vert	*point[3], *basevert;
	
	/* set source image globals */
	gfsource = fsource;
	gftarget = ftarget;
	gM1 = M1;
	gM2 = M2;

	for (i = 0; i < N2; i++) {
		gtargeti = i;
		for (j = 0; j < N1; j++) {
			gtargetj = j;
			gfsum = gareasum = 0.0;
			/* do the upper triangle */
			point[0] = makevertex(j, i + 1);
			point[1] = makevertex(j + 1, i + 1);
			point[2] = makevertex(j + 1, i);
			goodtriangle = 1;
			for (index = 0; index < 3; index++) {
				goodtriangle *= deflection(point[index]->y, point[index]->x, &di, &dj);
				point[index]->y += di;
				point[index]->x += dj;
			}
			if (goodtriangle) {
				makering(&basevert, point, 3);
				if (globalmapmode == INVERSEMAPMODE) {
					garea = trianglearea(basevert);
				}
				decompose(basevert);
			}
			/* do the lower triangle */
			point[0] = makevertex(j, i + 1);
			point[1] = makevertex(j, i);
			point[2] = makevertex(j + 1, i);
			goodtriangle = 1;
			for (index = 0; index < 3; index++) {
				goodtriangle *= deflection(point[index]->y, point[index]->x, &di, &dj);
				point[index]->y += di;
				point[index]->x += dj;
			}
			if (goodtriangle) {
				makering(&basevert, point, 3);
				if (globalmapmode == INVERSEMAPMODE) {
					garea = trianglearea(basevert);
				}
				decompose(basevert);
			}
			/* now add the pixel value */
			if ((gareasum > 0.0) && (globalmapmode == FORWARDMAPMODE))
				ftarget[i][j] = gfsum / gareasum;
			freeverts();
		}
	}
}



void	addarea(int ix, int iy, float area)
{
	if (ix >= 0 && ix < gM1 && iy >= 0 && iy < gM2) {
		if (globalmapmode == FORWARDMAPMODE) {
			if (gfsource[iy][ix] != FLOAT_MAGIC) {
				gareasum += area;
				gfsum += area * gfsource[iy][ix];
			}
		} else {
			if ((gftarget[gtargeti][gtargetj] != FLOAT_MAGIC) && (garea >= 0.0)) {
				gfsource[iy][ix] += area * gftarget[gtargeti][gtargetj] / garea;
			}
		}
	}
	gntriangles++;
}


long	trianglecount(void)
{
	return (gntriangles);
}



/* globals used by allocvert() and freeverts() */

static	int	counter = 0;
#define	MAX_VERTS	100000
static	vert	vstore[MAX_VERTS];

vert	*allocvert(void)
{
	counter++;
	if (counter >= MAX_VERTS)
		error_exit("allocvert: I'm out of verts\n");
	return (vstore + counter);
}


void	freeverts(void)
{
	counter = 0;
}



void	set_triangle_map_mode(int mapmode)
{
	globalmapmode = mapmode;
}

#undef TINY
#undef	NOCUT
#undef	XCUT
#undef	YCUT
#undef	MAX_VERTS
