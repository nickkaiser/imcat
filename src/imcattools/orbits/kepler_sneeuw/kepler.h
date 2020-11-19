/*
 * kepler.c - routines to convert between Kepler elements and cartesian r[], v[]
 * 
 * method due to Kaula 66 - from notes on www by Nico Sneeuw
 *
 * Nick Kaiser 11/25/02
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
void	keplertocartesian(keplerorbit *theorbit, double *r, double *v);
void	cartesiantokepler(double *r, double *v, keplerorbit *theorbit);
