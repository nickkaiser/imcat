/*
 * diff2Dpolymodel.c
 */


#define usage "\n\
NAME\n\
	diff2Dpolymodel --- differentiate a 2-dimensional spatial polynomial model\n\
\n\
SYNOPSIS\n\
	diff2Dpolymodel\n\
\n\
DESCRIPTION\n\
	'diff2Dpolymodel' reads from stdin a 2D polynomial model file with\n\
	parameters for a single variable and sends to stdout a 2D poly model\n\
	file for the two spatial derivatives of the variable.\n\
\n\
	This program requires that for each l, the mode values are\n\
	supplied in the standard order: i.e. 0,1, ..... l.\n\
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
	int 	nmodes_in, *l_in, *m_in, asize_in, nvar_in;
	int 	nm, m, *l_out, *m_out, asize_out;
	double	**a_in, **a_out;
	char	*vardef[MODEFUNC_MAX_VARS], *xname, fname[128];


	if (argc > 1) {
		error_exit(usage);
	}

	get2Dpolymodel(NULL, &l_in, &m_in, &asize_in, &a_in, &nmodes_in, &nvar_in, vardef, &xname);
	if (nvar_in != 3) {
		error_exit("diff2Dpolymodel: sorry, I can only handle single variable par files\n");
	}
	if (strncmp(vardef[2], "1 1", 3)) {
		error_exit("diff2Dpolymodel: sorry, I can only handle single variable par files\n");
	}
	sscanf(vardef[2] + 3, "%s", fname);
	sprintf(vardef[2], "1 2 %s", fname);
	asize_out = 2;
	a_out = (double **) calloc(2, sizeof(double *));
	a_out[0] = (double *) calloc(nmodes_in, sizeof(double));
	a_out[1] = (double *) calloc(nmodes_in, sizeof(double));
	l_out = (int *) calloc(nmodes_in, sizeof(int));
	m_out = (int *) calloc(nmodes_in, sizeof(int));

	nm = 0;
	for (m = 0; m < nmodes_in; m++) {
		if ((l_in[m] > 0) && (m_in[m] < l_in[m])) {
			l_out[nm] = l_in[m] - 1;
			m_out[nm] = m_in[m];
			a_out[0][nm] = (l_in[m] - m_in[m]) * a_in[0][m];
			a_out[1][nm] = m_in[m+1] * a_in[0][m+1];
			nm++;
		}
		
	}
	write2Dpolymodel(NULL, nm, l_out, m_out, asize_out, a_out, 1, vardef + 2, xname);

}

