#define usage "\n\n\n\
NAME\n\
	fit1object - fit FITS file to generalised Gaussian model\n\
\n\
SYNOPSIS\n\
	fit1object [options....]\n\
		-u		# print this message\n\
\n\
DESCRIPTION\n\
	\"fit1object\" reads a fits image f from stdin and fits\n\
	this to a simple model object of the form\n\
\n\
		f_model = f0 exp(-0.5 * [q_ij (r-d)_i (r-d)_j]^alpha)\n\
\n\
	by minimising\n\
\n\
		sum_r fabs(f - f_model)^beta\n\
\n\
	By default, alpha = 1 and beta = 2, so the program is\n\
	fitting a 2-dimensional gaussian by least squares.\n\
	and the ouput is in the form:\n\
		f0, d_x, d_y, q_xx, q_xy, q_yy\n\
	the matrix q_ij being taken to be symmetric.\n\
	You can modify the behaviour with the following parameters:\n\
		-a alpha	# set slope of exponential arg\n\
		-b beta		# type of fit\n\
		-i		# output inverse of q_ij\n\
		-A		# output a, b, phi (TBI)\n\
		-f		# output the model as fits image\n\
        The position is measured relative to the bottom left\n\
        corner of the bottom left pixel (so e.g. a single 'hot' pixel\n\
        at (ix,iy) = (23,67), would generate an object with x = (23.5, 67.5)\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@ifa.hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../../imlib/fits.h"
#include "../../../utils/error.h"
#include "gengaussfit.h"

#ifndef PI
#define PI M_PI
#endif

int		main(int argc, char *argv[])	
{
	int		arg = 1, N1, N2;
	int		invertq, outputabphi, outputmodelimage;
	float		**f;
	float		f0, d[2], **q;
	float		alpha, beta;
	float		det, temp;
/*	
	float		x, y, a, b, phi, f0;
*/
	fitsheader	*fits;

	/* defaults */
	alpha = 1.0;
	beta = 2.0;
	outputmodelimage = 0;
	invertq = 0;
	outputabphi = 0;

	/* parse args */	
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
			case 'a':
				sscanf(argv[arg++], "%f", &alpha);
				break;
			case 'b':
				sscanf(argv[arg++], "%f", &beta);
				break;
			case 'i':
				invertq = 1;
				break;
			case 'f':
				outputmodelimage = 1;
				break;
			case 'A':
				error_exit("fit1object: -A option not yet implemented\n");
				outputabphi = 1;
				invertq = 1;
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}

	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);

	/* initialise parameters */
	d[0] = 0.5 * N1;
	d[1] = 0.5 * N2;
	f0 = 0.0;
	q = (float **) calloc(2, sizeof(float *));
	q[0] = (float *) calloc(2, sizeof(float));
	q[1] = (float *) calloc(2, sizeof(float));

	/* do the fit */
	gaussfit(f, N1, N2, d, q, &f0, alpha, beta);

	/* output results */
	if (outputmodelimage) {
		add_comment(argc, argv, fits);
		makemodel(f, N1, N2, d, q, f0);
		write2Dfloatimage(f, fits);
		exit(0);
	}

	if (invertq) {
		det = (q[0][0] * q[1][1] - q[0][1] * q[1][0]);
		if (det == 0.0)
			error_exit("fit1object: zero determinent!\n");
		temp = q[0][0];
		q[0][0] = q[1][1] / det;
		q[1][1] = temp / det;
		q[0][1] = q[1][0] = - q[0][1] / det; 
	}

	fprintf(stdout, "%13g %13g %13g %13g %13g %13g\n", 
		f0, 0.5 + d[0], 0.5 + d[1], q[0][0], q[0][1], q[1][1]);
	exit(0);
}



