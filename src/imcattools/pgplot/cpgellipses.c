/*
 * pgellipse.c --- code to draw ellipses with pgplot
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cpgplot.h"


#define PI M_PI

/*
 cpgellipses() -- draw a set of ne ellipses at locations x, y with major,minor axes a, b and pos angle phi0 (rad)
	using a polygon with nv vertices per quadrant, so nv = 1 gives a diamond, nv = 2 has 8 sides etc...
	If colour != NULL we fill the polygon first with a colour.
	If colour == NULL and fill >= 0.0 we fill polygon with shade = fill.
	If drawoutline == 1 we draw the outline of the ellipse last in foreground colour
*/
void	cpgellipses(int ne, float *x, float *y, float *a, float *b, float *phi0, float *colour, 
		int nv, int drawoutline, float fill)
{
	static	int	lastnv = 0;
	static	float	*xv = NULL, *yv = NULL;
	float	X, Y, c, s, phi, dphi;
	int	e, v, cihi, cilo, oldci, newci;

	if (!nv) {
		fprintf(stderr, "pgellipse: non-positive nv!\n");
		exit(-1);
	}
	if (nv != lastnv) {		/* we need to allocate the vertices */
		if (xv)
			free(xv);
		if (yv)
			free(yv);
		xv = (float *) calloc(4 * nv + 1, sizeof(float));
		yv = (float *) calloc(4 * nv + 1, sizeof(float));
		lastnv = nv;
	}

	cpgqcir(&cilo, &cihi);
	
	dphi = 0.5 * PI / nv;
	for (e = 0; e < ne; e++) {		/* loop over ellipses */
		c = cos(phi0[e]);
		s = sin(phi0[e]);
		phi = 0.0;
		/* make the ellipse */
		for (v = 0; v < 4 * nv + 1; v++) {
			X = a[e] * cos(phi);
			Y = b[e] * sin(phi);
			xv[v] = x[e] + X * c - Y * s;
			yv[v] = y[e] + X * s + Y * c;
			phi += dphi;
		}
		if (colour || (fill >= 0.0)) {
			cpgqci(&oldci);
			if (colour) {
				newci = cilo + (int) floor(colour[e] * (1 + cihi - cilo));
			} else {
				newci = cilo + (int) floor(fill * (1 + cihi - cilo));
			}
			newci = (newci > cihi ? cihi : (newci < cilo ? cilo : newci));
			cpgsci(newci);
			cpgpoly(4 * nv, xv, yv);
			cpgsci(oldci);
		}
		if (drawoutline) {
			cpgline(4 * nv + 1, xv, yv);
		}
	}			
}
