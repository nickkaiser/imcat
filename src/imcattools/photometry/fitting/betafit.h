/*
 * betafit.h
 */

float	func(float *p);
int     fitall(float **n, int N1, int N2, float *f0, float *x0, float *y0, float *rc2, float beta, float *loglhood, float N);
void	makebetamodel(float **nmodel, int N1, int N2, float f0, float x0, float y0, float rc2, float beta, float N);
