/*
 * kepler elements to cartesian conversion using Heasley's routines
 */

typedef struct keplerorbit {
	double	a;
	double	e;
	double	i;
	double	omega;
	double	Omega;
	double	M;
} keplerorbit;


void	assignkeplerorbitfrombuffer(keplerorbit *orbit, double *buffer);
void 	fillbufferfromkeplerorbit(keplerorbit *orbit, double *buffer);
int	keplertocartesian(keplerorbit *theorbit, double *r, double *v);
int	cartesiantokepler(double *r, double *v, keplerorbit *theorbit);
