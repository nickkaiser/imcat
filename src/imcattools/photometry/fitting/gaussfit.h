/*
 * gaussfit.h
 */

float	func(float *p);
void	dfunc(float *p, float *df);
int     gaussfit(float **ff, int N1, int N2, float *f0, float *d, float *a, float *b, float *phi);
void	makemodel(float **ff, int N1, int N2, float f0, float *d, float a, float b, float phi);

