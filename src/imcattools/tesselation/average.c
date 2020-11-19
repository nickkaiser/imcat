/*
 * average.c
 */

#include <math.h>


double	mean(double **y, int j)
{
	int 	c;
	double	res;

	res = 0;
	for (c = 0; c < 3; c++) {
		res += y[c][j];
	}
	return(res / 3.0);
}

double	median(double **y, int j)
{
	if (((y[1][j] - y[0][j]) * (y[2][j] - y[0][j])) < 0.0) {
		return(y[0][j]);
	}
	if (((y[2][j] - y[1][j]) * (y[0][j] - y[1][j])) < 0.0) {
		return(y[1][j]);
	}
	return(y[2][j]);
}

