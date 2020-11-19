/*
 * difflmodel.c
 */


#define usage "\n\
NAME\n\
	difflmodel --- differentiate linear mode superposition model\n\
\n\
SYNOPSIS\n\
	difflmodel [-i]\n\
\n\
DESCRIPTION\n\
	'difflmodel' reads from stdin a catalogue containing the\n\
	definition of a linear mode function superposition model\n\
	or 'lmodel' and sends to stdout a lmodel giving the\n\
	coefficients of the derivative of the input model function.\n\
\n\
	For example, given a model of a rank-2 matrix\n\
	valued function f_ij(x) on an n-dimensional space\n\
	which we model as:\n\
\n\
		F_ij(x) = sum_m a_mij f_m(x)\n\
\n\
	the result is a set of mode coefficients a'_mijlm such that\n\
\n\
		F'_ijl(x) = sum_m a_mijl f_m(x) = d F_ij / d x_l\n\
\n\
	With -i option we perform the inverse operation - i.e. we\n\
	compute the integral of the input model.\n\
\n\
	It only works for polynomial or Fourier models.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "catlib/cat.h"
#include "utils/error.h"
#include "utils/args.h"
#include "utils/lmodel.h"


main(int argc, char *argv[])
{
	lmodel	*ipmodel, *opmodel;
	char	*flag;
	int	inversemode = 0;

	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'i':
				inversemode = 1;
				break;
			default:
				error_exit(usage);
		}
	}

	/* read the source lmodel */
	ipmodel = readlmodel(stdin);

	/* do the biz */
	opmodel = difflmodel(ipmodel, inversemode);

	/* output the result */
	writelmodel(opmodel, stdout);

	exit(0);
}


