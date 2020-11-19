/*
 * kernels.c
 *
 * functions to calculate kernels for simple disc or rectangle case b
 * and cutoff(a,eps) which does some cut-off
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernels.h"

#define EPS 1.e-10


/* next function calculates K = kernel for a reference disk of radius R */
/* the target point r0 and galaxy position r can lie outside the disk */
/* eps is the cut-off parameter to stop things blowing up */
void	circkernel(double *K, double x0, double y0, double x, double y, double R, double eps)
{
	double	x2, y2, r, r2, r0, r02, R2, mu, z, zmu, Y, KK, a1, a2, W;

	/* convert abs to rel measurement pos */
	x -= x0;
	y -= y0;

	/* calculate mu, z = discriminant */
	x2 = x * x;
	y2 = y * y;
	r2 = x2 + y2 + EPS;
	r = sqrt(r2);
	r02 = x0 * x0 + y0 * y0 + EPS;
	r0 = sqrt(r02);
	R2 = R * R;
	mu = (x * x0 + y * y0) / (r * r0);
	z = R2 / r02 - 1 + mu * mu;

	/* check to see if we don't intersect disk at all: if so set K = 0 */
	if (z <= 0.0) {
		K[0] = K[1] = 0.0;
		return;
	}

	/* discriminant is positive so we can solve for a1 = p1 / r, a2 = p2 / r,  a2 > a1 */
	z = sqrt(z);
	a1 = - r0 * (mu + z) / r;
	a2 = - r0 * (mu - z) / r;

	/* check to see if we're in front of disk and r is behind: if so set K = 0 */
	if (a2 <= 1.0) {
		K[0] = K[1] = 0.0;
		return;
	}

	/* deal with the largest root a2 */
	zmu = z - mu;
	W = cutoff(r, a2 * r, eps);
	KK = - W * r02 * zmu * zmu / (z * r2 * r2);
	K[0] = KK * (zmu * (x2 - y2) + r * (x * x0 - y * y0) / r0);
	K[1] = KK * (zmu * 2 * x * y + r * (x * y0 + y * x0) / r0);

	/* check to see if 2nd root is +ve; if so target point is outside disk so subtract */
	if (a1 > 0) {
		zmu = z + mu;
		W = cutoff(r, a1 * r, eps);
		KK = - W * r02 * zmu * zmu / (z * r2 * r2);
		K[0] -= KK * (zmu * (x2 - y2) - r * (x * x0 - y * y0) / r0);
		K[1] -= KK * (zmu * 2 * x * y - r * (x * y0 + y * x0) / r0);
	}
}


#define	sign(x)	((x) > 0 ? 1 : -1)

/* calculate the kernel for a rectangular reference region */
/* target point = (x0,y0), measurement point = (x,y) reference region is a */
/* box of half-side lengths LX, LY */
/* kernel = 0 for r < eps * p */
/* also returns alpha = r / p (set to zero for x0,y0 outside box) */
void	rectkernel(double *K, double *alpha, double x0, double y0, double x, double y, double LX, double LY, double eps)
{
	double	P, LL[2], L[2], rr0[2], r0[2], r[2], rr[2], a[2][2], aa[2], KK[2][2], W, modr;
	int	i, j, count;

	rr0[1] = r0[0] = x0;
	rr0[0] = r0[1] = y0;
	rr[1] = r[0] = x - x0;
	rr[0] = r[1] = y - y0;
	LL[1] = L[0] = LX;
	LL[0] = L[1] = LY + EPS;

	/* get four solutions for signed a^2 = p^2 / r^2 */
	/* these give intersections of the line with the lines x = +- LX, y = +- LY */
	for (i = 0; i < 2; i++) {
		a[i][0] = - (r0[i] + L[i]) * fabs(r0[i] + L[i]) * sign(r[i]) / (r[i] * r[i] + EPS);
		a[i][1] = - (r0[i] - L[i]) * fabs(r0[i] - L[i]) * sign(r[i]) / (r[i] * r[i] + EPS);
	}

	count = 0;
	K[0] = K[1] = 0.0;
	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			/* reject negative a solutions */
			if (a[i][j] <= 0.0)
				continue;
			/* calculate intersection with perimeter */
			P = rr0[i] + sqrt(a[i][j]) * rr[i];
			/* reject intersections outside the box */
			if (fabs(P) >= LL[i])
				continue;
			/* reject points with p < r */
			if (a[i][j] <= 1.0)
				continue;
			if (!i) {
				KK[count][0] = - a[i][j];
				KK[count][1] = r[0] * r[1] * KK[count][0] / (r[0] * r[0] + EPS);	
			} else {
				KK[count][0] = a[i][j];	
				KK[count][1] = - r[0] * r[1] * KK[count][0] / (r[1] * r[1] + EPS);	
			}
			/* multiply by  smoothing factor */
			modr = sqrt(r[0] * r[0] + r[1] * r[1]);
			W = cutoff(modr, modr * sqrt(a[i][j]), eps);
			KK[count][0] *= W;
			KK[count][1] *= W;
			aa[count] = a[i][j];
			count++;
			if (count == 2)		/* deals with special cases */
				break;
		}
		if (count == 2)
			break;
	}
	if (count) {
		K[0] = KK[0][0];
		K[1] = KK[0][1];
		*alpha = 1.0 / sqrt(aa[0]);
	}
	if (count > 1) {
		if (aa[1] > aa[0]) {
			K[0] -= KK[1][0];
			K[1] -= KK[1][1];
		} else {
			K[0] = KK[1][0] - K[0];
			K[1] = KK[1][1] - K[1];
		}
		*alpha = 0.0;
	}
}



/*
 * cutoff() implements small scale smoothing
 *
 * currently a function only of r / eps * p, but give args separately so we
 * can put in a cut off at certain number of pixels if necessary
 *
 */
double	cutoff(double r, double p, double eps)
{
	double	x;

	if (p < EPS)
		p = EPS;
	x = r / (eps * p);
	return((x < 2.0 ? 1.0 - exp(- x * x * x) : 1.0));
}


#undef	EPS

