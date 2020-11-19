/**
 ** function for testing Gaussian ellipsoid fitting
 **/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "gaussfit.h"

#ifndef PI
#define PI M_PI
#endif

main()
{
	float x, y, a, b, phi, f0;
	int	N = 10, nfit, i, j, iO, jO;
	float **f;
	float	Cx, Cy, Cxx, Cxy, Cyy, c, s, Q;
	
	f = (float **) calloc(N, sizeof(float *));
	for (i = 0; i < N; i++)
		f[i] = (float *) calloc(N, sizeof(float));
	printf("f0 x y?\n");
	scanf("%f %f %f", &f0, &x, &y);
	iO = floor(0.5 + x);
	jO = floor(0.5 + y);
	printf("a b phi?\n");
	scanf("%f %f %f", &a, &b, &phi);
	phi *= 2 * PI / 360.0;
	Cx = pow((double) a, -2.0);
	Cy = pow((double) b, -2.0);
	c = cos((double) phi);
	s = sin((double) phi);
	Cxx = Cx * c * c + Cy * s * s;
	Cyy = Cx * s * s + Cy * c * c;
	Cxy = (Cx - Cy) * c * s;
	
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			Q = Cxx*(i-x)*(i-x)+2*Cxy*(i-x)*(j-y)+Cyy*(j-y)*(j-y); 
			f[i][j] = f0 * exp(-0.5 * Q);
		}	
	}
	nfit = 3;
	
	gaussfit(f, nfit, iO, jO, &x, &y, &a, &b, &phi, &f0);
	fprintf(stdout, "\nx=%f\ny=%f\nf0=%f\n", iO + x, jO + y, f0);
	fprintf(stdout, "a=%f\nb=%f\nphi=%f\n\n", a, b, 360 * phi / (2 * PI));
}




