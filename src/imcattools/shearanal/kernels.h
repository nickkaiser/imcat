/*
 * kernels.h
 *
 * declaration of kernels for disk and rectangle - case-b only
 * call twice with different boundaries and subtract for case-c
 *
 * cutoff() implements smoothing
 *
 */

void	circkernel(double *K, double x0, double y0, double x, double y, double R, double eps);
void    rectkernel(double *K, double *alpha, double x0, double y0, double x, double y, double NX, double NY, double eps);
double	cutoff(double r, double p, double eps);


