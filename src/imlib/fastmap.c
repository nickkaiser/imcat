/*
 * fastmap.c - defines fastmap()
 */

#include <math.h>
#include <stdio.h>
#include <limits.h>
#include "map.h"
#include "fits.h"

#define MAGIC FLOAT_MAGIC

#define	W	0
#define	S	1
#define N	2
#define E	3

static int ix[8] = {0,0,0,1,0,1,1,1};
static int iy[8] = {0,1,0,0,1,1,0,1}; 
 
/*
 * fastmap() maps a M2 x M1 source image fsource[i][j] onto an N2 x N1 target 
 * image ftarget[i][j] with the mapping 
 *		ftarget(r) += fsource(r + d)
 * where the deflection d = (di, dj) is supplied by function deflection()
 * this version simply finds the nearest pixel
 */
void	ultrafastmap(float **ftarget, int N1, int N2, float **fsource, int M1, int M2,
			int (*deflection)(float ri, float rj, float *di, float *dj))
{
	int	i, j, is, js;
	float	x, y, f0, fx, fy, xs, ys;

	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {
			/* map the centre of the target pixel */
			if (deflection(i + 0.5, j + 0.5, &y, &x)) {
				xs = j + 0.5 + x;
				ys = i + 0.5 + y;
				/* get the pixel whose centre it lies closest to */
				is = (int) floor(ys);
				js = (int) floor(xs);
				if (is < 0 || is >= M2 || js < 0 || js >= M1)
					continue;
				ftarget[i][j] = fsource[is][js];
			}
		}
	}
}

/*
 * this version finds the quadrant (NSEW) that the point lives in,
 * and then does a linear interpolation using plane defined by the
 * centre of the source pixel and the two corners of the quadrant
 */


void	fastmap(float **ftarget, int N1, int N2, float **f, int M1, int M2,
			int (*deflection)(float ri, float rj, float *di, float *dj))
{
	int	i, j, is, js, q, ind;
	float	x, y, xs, ys, dx, dy, f0, fx, fy, f1, f2, fxg, fyg;

	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {
			/* map the centre of the target pixel */
			deflection(i + 0.5, j + 0.5, &y, &x);
			/* get the pixel it lives in */
			xs = j + x;
			ys = i + y;
			is = (int) floor(ys);
			js = (int) floor(xs);
			if (is < 0 || is >= (M2-1) || js < 0 || js >= (M1-1))
				continue;
			if ((f[is][js] == MAGIC) || (f[is][js+1] == MAGIC) || (f[is+1][js] == MAGIC) || (f[is+1][js+1] == MAGIC)) {
				ftarget[i][j] = MAGIC;
			} else {
				/* calculate the mean of the corner values */
				f0 = 0.25 * (f[is][js] + f[is][js+1] + f[is+1][js] + f[is+1][js+1]);
				/* and the 'global' gradients */
				fxg = 0.5 * (f[is][js+1] - f[is][js] + f[is+1][js+1] - f[is+1][js]);
				fyg = 0.5 * (f[is+1][js] - f[is][js] + f[is+1][js+1] - f[is][js+1]);
				/* calculate the displacement relative to the centre of source pixel */
				dx = xs - js - 0.5;
				dy = ys - is - 0.5;
				/* get the quadrant */
				q = 2 * ((dx + dy) > 0.0) + ((dx - dy) > 0);
				/* get the vertices of the quadrant */
				ind = 2 * q;
				f1 = f[is + iy[ind]][js + ix[ind++]];
				f2 = f[is + iy[ind]][js + ix[ind]];
				/* calculate the relevant gradients */
				switch (q) {
					case N:
					case S:
						fx = f2 - f1;
						fy = fyg;
						break;
					case E:
					case W:
						fy = f2 - f1;
						fx = fxg;
						break;
					default:
						fprintf(stderr, "fastmap: bad quadrant\n");
						exit(-1);
						break;
				}	
				ftarget[i][j] = f0 + dx * fx + dy * fy;
			}
		}
	}
}

