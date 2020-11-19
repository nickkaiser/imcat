/*
 * colormaps.c
 */

#include <stdlib.h>
#include <stdio.h>

#define CMAP_INDEX_MAX 5


int	getcolormap(float **lp, float **rp, float **gp, float **bp, float *contra, float *bright, int *n, int index)
{
	float	*l, *r, *g, *b;
	int	c;
	FILE	*cmapf;
	char	filename[512];
	char	*cmap[CMAP_INDEX_MAX];

	cmap[0] = "greece.cmap";
	cmap[1] = "candy.cmap";
	cmap[2] = "pretty.cmap";
	cmap[3] = "terra.cmap";
	cmap[4] = "ocean.cmap";

	if (index >= CMAP_INDEX_MAX) {
		error_exit("getcolormap: illegal index\n");
	}
	sprintf(filename, "%s", getenv("IMCATDIR"));
	strcat(filename, "/src/utils/colormaps/");
	strcat(filename, cmap[index]);
	cmapf = fopen(filename, "r");
	if (!cmapf) {
		error_exit("getcolormap: failed to open colormap\n");
	}
	fscanf(cmapf, "%d", n);
	l = (float *) calloc(*n, sizeof(float));
	r = (float *) calloc(*n, sizeof(float));
	g = (float *) calloc(*n, sizeof(float));
	b = (float *) calloc(*n, sizeof(float));
	for (c = 0; c < *n; c++) {
		l[c] = ((float) c) / (*n - 1);
		fscanf(cmapf, "%f %f %f", r + c, g + c, b + c);
	}
	*lp = l;
	*rp = r;
	*gp = g;
	*bp = b;
	*contra = 1.0;
	*bright = 0.5;
}

int	getrgbfromcmap(float **r, float **g, float **b, int n, int cmapindex)
{
	int 	i, j, m;
	float 	contra, bright, *l, *rr, *gg, *bb, f, frac;

	getcolormap(&l, &rr, &gg, &bb, &contra, &bright, &m, cmapindex);
	*r = (float *) calloc(n, sizeof(float));
	*g = (float *) calloc(n, sizeof(float));
	*b = (float *) calloc(n, sizeof(float));

	for (i = 0; i < n; i++) {
		f = (n - i - 0.5) / n;
		for (j = 0; j < m; j++) {
			if (l[j] > f) {
				break;
			}
		}
		j = (j == m ? m - 1 : j);
		j = (j == 0 ? 1 : j);
		frac = (l[j] - f) / (l[j] - l[j - 1]);
		(*r)[i] = (1 - frac) * rr[j] + frac * rr[j - 1];
		(*g)[i] = (1 - frac) * gg[j] + frac * gg[j - 1];
		(*b)[i] = (1 - frac) * bb[j] + frac * bb[j - 1];
	}
}




