/*
 * average.h
 */

/* f-averaging modes */
#define F_INTERP	0
#define	F_MEAN		1
#define F_MEDIAN	2

double	mean(double **y, int j);
double	median(double **y, int j);
