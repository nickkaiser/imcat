/*
 * gaussfit.h
 */

float	func(float *p);
void	dfunc(float *p, float *df);
int     gaussfit(float **ff, int N1, int N2, float *d, float **q, float *f0,
        float alpha, float beta);
void	makemodel(float **ff, int N1, int N2, float *d, float **q, float f0);

