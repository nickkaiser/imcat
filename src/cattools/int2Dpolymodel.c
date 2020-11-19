/*
 * int2Dpolymodel.c
 */


#define usage "\n\
NAME\n\
	int2Dpolymodel --- integrate a 2-dimensional spatial polynomial model\n\
\n\
SYNOPSIS\n\
	int2Dpolymodel\n\
\n\
DESCRIPTION\n\
	'int2Dpolymodel' reads from stdin a 2D polynomial model file with\n\
	parameters for a 2-vector variable (which is supposed to represent\n\
	the spatial derivative of some scalar function 'f' [thought the\n\
	actual name is arbitrary]) and sends to stdout a 2D poly model\n\
	file for the mode coefficients of f. More explicitly, given a\n\
	set of pairs of mode amplitudes a^0_{l,m} a^1_{l,m}  then we output\n\
	a single set of mode amplitudes:\n\
		a_{l,m}  = ((l-m) a^0_{l-1,m} + m a^1_{l-1,m-1}) / (m^2 + (l-m)^2)\n\
	which is the optimal combination assuming statistically independent,\n\
	but equally noisy, coefficients for a^0, a^1.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../utils/modefunc.h"

main(int argc, char *argv[])
{
	int 	nmodes_in, *l_in, *m_in, asize_in, nvar_in, lmin, lmax, l0, l1, m0, m1;
	int 	nmodes_out, l, m, *l_out, *m_out, asize_out, **mode;
	double	**a_in, **a_out, lo, mo, denom;
	char	*vardef[MODEFUNC_MAX_VARS], *xname, fname[128];
	int	nm;


	if (argc > 1) {
		error_exit(usage);
	}

	get2Dpolymodel(NULL, &l_in, &m_in, &asize_in, &a_in, &nmodes_in, &nvar_in, vardef, &xname);
	if (nvar_in != 3) {
		error_exit("int2Dpolymodel: sorry, I can only handle single variable par files\n");
	}
	if (strncmp(vardef[2], "1 2", 3)) {
		error_exit("int2Dpolymodel: mode amplitudes must be a 2 vector\n");
	}

	/* we need to figure lmin, lmax and create mode[l][m] */
	lmax = 0;
	lmin = 100000;
	for (m = 0; m < nmodes_in; m++) {
		if (l_in[m] > lmax) {
			lmax = l_in[m];
		}
		if (l_in[m] < lmin) {
			lmin = l_in[m];
		}
	}
	mode = (int **) calloc(lmax + 1, sizeof(int *));
	for (l = lmin; l <= lmax; l++) {
		mode[l] = (int *) calloc(lmax + 1, sizeof(int));
		for (m = 0; m <= lmax; m++) {
			mode[l][m] = -1;
		}
	}
	for (m = 0; m < nmodes_in; m++) {
		mode[l_in[m]][m_in[m]] = m;
	}

	/* compute nmodes_out */
	lmin++;
	lmax++;
	nmodes_out = 0;
	for (l = lmin; l <= lmax; l++) {
		nmodes_out += (l + 1);
	}
	/* allocate l_out, m_out, a_out */
	a_out = (double **) calloc(1, sizeof(double *));
	a_out[0] = (double *) calloc(nmodes_out, sizeof(double));
	l_out = (int *) calloc(nmodes_out, sizeof(int));
	m_out = (int *) calloc(nmodes_out, sizeof(int));
	/* compute l_out, m_out */
	nm = 0;
	for (l = lmin; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			l_out[nm] = l;
			m_out[nm] = m;
			nm++;
		}
	}

	/* now we do the biz */
	for (nm = 0; nm < nmodes_out; nm++) {
		lo = l_out[nm];
		mo = m_out[nm];
		denom = mo * mo + (lo - mo) * (lo - mo);
		l0 = l1 = l_out[nm] - 1;
		if (m_out[nm] < l_out[nm]) {
			m0 = m_out[nm];
			a_out[0][nm] += (lo - mo) * a_in[0][mode[l0][m0]] / denom;
		}
		if (m_out[nm] > 0) {
			m1 = m_out[nm] - 1;
			a_out[0][nm] += mo * a_in[1][mode[l1][m1]] / denom;
		}
	}

	sscanf(vardef[2] + 3, "%s", fname);
	sprintf(vardef[2], "1 1 %s", fname);
	asize_out = 1;
	write2Dpolymodel(NULL, nm, l_out, m_out, asize_out, a_out, 1, vardef + 2, xname);

}

