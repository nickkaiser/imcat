#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../imlib/fits.h"

unsigned char ***makecarray(fitsheader *fits, int N1, int N2, int Ncolors, float fmin, float fmax, int dobox) 
{
	float		*f, fval;
	unsigned char 	***c;
	int		color, i, j;

	f = (float *) calloc(N1, sizeof(float));
	c = (unsigned char ***) calloc(Ncolors, sizeof(unsigned char **));
	for (color = 0; color < Ncolors; color++) {
		c[color] = (unsigned char **) calloc(N2, sizeof(unsigned char *));
		for (i = 0; i < N2; i++) {
			c[color][i] = (unsigned char *) calloc(N1, sizeof(unsigned char));
			readfitsline(f, fits);
			for (j = 0; j < N1; j++) {
				if (Ncolors != 1) {
					if (f[j] == FLOAT_MAGIC) {
						fval = 255;
					} else {
						fval = 256 * (f[j] - fmin) / (fmax - fmin);
					}
				} else {
					fval = 256 * (fmax - f[j]) / (fmax - fmin);
				}
				c[color][i][j] = (unsigned char) (fval >= 0 ? (fval < 256 ? fval : 255) : 0);
			}
		}
		if (dobox) {
			for (i = 0; i < N1; i++) {
				c[color][0][i] = c[color][N2 - 1][i] = (unsigned char) 0;
			}
			for (i = 0; i < N2; i++) {
				c[color][i][0] = c[color][i][N1 - 1] = (unsigned char) 0;
			}
		}
	}
	return(c);
}
