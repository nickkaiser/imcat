/*
 * colormaps.h
 */

int	getcolormap(float **lp, float **rp, float **gp, float **bp, float *contra, float *bright, int *n, int index);

int getrgbfromcmap(float **r, float **g, float **b, int n, int cmapindex);

