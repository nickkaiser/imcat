/*
 * efit_stuff.h
 */

/*
 * declaration for mode function g(mode, frame, xi, xj) 
 * and function for reading and writing mode amplitudes
 *
 * the model ellipticity is given by
 *	e[pol] = sum_modes a[pol][mode] * g(mode, frame, i, j)
 * where i, j are pixel coords
 *
 * don't forget to call setnframes() and setframsize() before calling g()
 */

#define	GRAD1 0
#define	GRAD2 1
#define	GRADGRAD11 2
#define	GRADGRAD12 3
#define	GRADGRAD22 4

double		g(int mode, int frame, double xx, double yy);
void		readamplitudes(int *nmodes, int *nframes, int *framesize, 
			int *order, double ***amp, FILE *stream);
void		setnframes(int nframes);
void		setframesize(int N);

