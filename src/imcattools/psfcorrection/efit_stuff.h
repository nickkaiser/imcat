/*
 * efit_stuff.h
 */

/*
 * declaration for mode function g(mode, xi, xj) 
 * and functions for reading and writing mode amplitudes
 *
 * the model ellipticity is given by
 *	e[pol] = sum_modes a[pol][mode] * g(mode, frame, x, y)
 * where x, y are coords
 *
 * don't forget to call setframsize() before calling g()
 */


double		g(int mode, double xx, double yy);
void		readamplitudes(int *nmodes, int *framesize, int *order, double ***amp, FILE *stream);
void		writeamplitudes(int nmodes, int order, int framesize, double **amp);
void		setframesize(int N);

