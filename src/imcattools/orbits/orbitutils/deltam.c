#include <stdio.h>
#include <math.h>
#include "deltam.h"

double  deltam(double *rhelio, double *rgeo, double G)
{
	double	A[2], B[2], rr, dd, dr, tiny, alpha;
	double	phi[2], dm;
	int	l;

	A[0] = 3.33; 
	A[1] = 1.87;
	B[0] = 0.63;
	B[1] = 1.22;
	tiny = 1.e-20;

	rr = sqrt(rhelio[0] * rhelio[0] + rhelio[1] * rhelio[1]);
	dd = sqrt(rgeo[0] * rgeo[0] + rgeo[1] * rgeo[1]);
	dr = rhelio[0] * rgeo[0] + rhelio[1] * rgeo[1];
	if ((rr > tiny) && (dd > tiny) && (fabs(dr) > tiny)) {
		alpha = acos(dr / (rr * dd));
		for (l = 0; l < 2; l++) {
			phi[l] = exp(- A[l] * pow(tan(alpha / 2.0), B[l]));
		}
		dm = - 2.5 * log10((1 - G) * phi[0] + G * phi[1]) + 5.0 * log10(dd * rr);
	} else {
		dm = 0.0;
	}

	return(dm);
}
