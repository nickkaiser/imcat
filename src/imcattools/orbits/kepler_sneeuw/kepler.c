/*
 * kepler.c - routines to convert between Kepler elements and cartesian r[], v[]
 * 
 * method due to Kaula 66 - from notes on www by Nico Sneeuw 
 * 
 * later replaced by www.bruce-shapiro.net/pair/ ElementConversionRecipes.pdf 
 *
 * Nick Kaiser 11/25/02
 */

#include <stdio.h>
#include <math.h>
#include "vectors.h"
#include "kepler.h"

#define KEPLER_TINY 1.e-15

void 	assignkeplerorbitfrombuffer(keplerorbit *orbit, double *buffer)
{
	orbit->a = buffer[0];
	orbit->e = buffer[1];
	orbit->i = M_PI * buffer[2] / 180.0;
	orbit->omega = M_PI * buffer[3] / 180.0;
	orbit->Omega = M_PI * buffer[4] / 180.0;
	orbit->M = M_PI * buffer[5] / 180.0;
}

void 	fillbufferfromkeplerorbit(keplerorbit *orbit, double *buffer)
{
	buffer[0] = orbit->a;
	buffer[1] = orbit->e;
	buffer[2] = 180.0 * orbit->i / M_PI;
	buffer[3] = 180.0 * orbit->omega / M_PI;
	buffer[4] = 180.0 * orbit->Omega / M_PI;
	buffer[5] = 180.0 * orbit->M / M_PI;
}

/* from Sneeuw (except for a^{-3/2} factor in v) */
void	keplertocartesian(keplerorbit *theorbit, double *r, double *v)
{
	double	E, Enew, a, e, i, omega, Omega, M;

	a = theorbit->a; 
	e = theorbit->e; 
	i = theorbit->i; 
	omega = theorbit->omega; 
	Omega = theorbit->Omega; 
	M = theorbit->M;

	/* solve iteratively for the eccentric anomaly */
	E = M;
	while (1) {
		Enew = e * sin(E) + M;
		if (fabs(Enew -  E) < 1.e-10) {
			break;
		}
		E = Enew;
	}

	assign(r, a * (cos(E) - e), a * sqrt(1 - e * e) * sin(E), 0.0);
	assign(v, - sin(E), sqrt(1 - e * e) * cos(E), 0.0);
	if ((1 - e * cos(E)) > KEPLER_TINY) {
		scale(v, 1.0 / (sqrt(a) * (1 - e * cos(E))));
	}
	rotz(r, - omega);
	rotx(r, - i);
	rotz(r, - Omega);
	rotz(v, - omega);
	rotx(v, - i);
	rotz(v, - Omega);

	/* this is a kludge to bring things into accord with Shapiro */
	r[1] *= -1.0;
	v[0] *= -1.0;
	v[2] *= -1.0;
}


void	keplertocartesian_shapiro(keplerorbit *theorbit, double *r, double *v)
{
	double	E, Enew, a, e, i, omega, Omega, M, *P, *Q, A, B;
	int	j;

	a = theorbit->a; 
	e = theorbit->e; 
	i = theorbit->i; 
	omega = theorbit->omega; 
	Omega = theorbit->Omega; 
	M = theorbit->M;

	/* solve iteratively for the eccentric anomaly */
	E = M;
	while (1) {
		Enew = e * sin(E) + M;
		if (fabs(Enew -  E) < KEPLER_TINY) {
			break;
		}
		E = Enew;
	}

	P = (double *) calloc(3, sizeof(double));
	Q = (double *) calloc(3, sizeof(double));
	assign(P, 
		cos(omega) * cos(Omega) - sin(omega) * cos(i) * sin(Omega),
		cos(omega) * sin(Omega) + sin(omega) * cos(i) * cos(Omega),
		sin(omega) * sin(i));
	assign(Q,
		-sin(omega) * cos(Omega) - cos(omega) * cos(i) * sin(Omega),
		-sin(omega) * sin(Omega) + cos(omega) * cos(i) * cos(Omega), 
		sin(i) * cos(omega));
	
	A = a * (cos(E) - e);
	B = a * sqrt(1 - e * e) * sin(E);
	for (j = 0; j < 3; j++) {
		r[j] = A * P[j] + B * Q[j];
	}
	A = - e * sin(E) / ((1 - e * cos(E)) * sqrt(a));
	B = sqrt(1 - e * e) * cos(E) / ((1 - e * cos(E)) * sqrt(a));
	for (j = 0; j < 3; j++) {
		v[j] = A * P[j] + B * Q[j];
	}
}

/* from Shapiro */
void	cartesiantokepler(double *r, double *v, keplerorbit *theorbit)
{
	double	a, i, omega, Omega, M;
	double	*h, modh, *e, mode, modv, modr, *rhat, *k, *n, modn, theta, E;
	int	j;

	h = (double *) calloc(3, sizeof(double));
	e = (double *) calloc(3, sizeof(double));
	rhat = (double *) calloc(3, sizeof(double));
	k = (double *) calloc(3, sizeof(double));
	n = (double *) calloc(3, sizeof(double));

	/* k[] is the unit z-vector */
	assign(k, 0.0, 0.0, 1.0);

	/* compute angular momentum h[] */
	cross(h, r, v);
	modh = length(h);
	
	/* compute the ellipticity vector e[] */
	cross(e, v, h);
	modr = length(r);
	copy(rhat, r);
	scale(rhat, 1.0 / modr);
	for (j = 0; j < 3; j++) {
		e[j] -= rhat[j];
	}
	mode = length(e);

	/* semi-major axis */
	a = modh * modh / (1 - mode * mode);

	/* inclination */
	i = acos(h[2] / modh);

	/* compute n = k x h */
	cross(n, k, h);
	modn = length(n);

	/* RA of ascending node Omega */
	Omega = acos(n[0] / modn);
	if (n[1] < 0.0) {
		Omega = 2.0 * M_PI - Omega;
	}

	/* argument of perigee */
	omega = acos(dot(n, e) / (modn * mode));
	if (dot(e, k) < 0.0) {
		omega = 2.0 * M_PI - omega;
	}

	/* true anomaly theta */
	theta = acos(dot(e, r) / (mode * modr));
	if (dot(r, v) > 0.0) {
		theta = 2 * M_PI - theta;
	}

	/* eccentric anomaly E */
	E = acos((mode + cos(theta)) / (1.0 + mode * cos(theta)));
	if ((M_PI < theta) && (theta < 2 * M_PI)) {
		E = 2 * M_PI - E;
	}

	/* mean anomaly */
	M = E - mode * sin(E);

	/* stack the keplerorbit structure */
	theorbit->a = a;
	theorbit->e = mode;
	theorbit->i = i;
	theorbit->omega = omega;
	theorbit->Omega = Omega;
	theorbit->M = M;

	free(h);
	free(e);
	free(rhat);
	free(k);
	free(n);

}

