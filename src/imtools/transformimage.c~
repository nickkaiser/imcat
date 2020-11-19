/*
 * transformimage.c
 */

#define	usage "\n\n\n\
NAME\n\
	transformimage - apply spatial linear transformation to a FITS image\n\
\n\
SYNOPSIS\n\
	transformimage [options...]\n\
		-p psi_xx psi_xy psi_yx psi_yy	# distortion matrix (1,0,0,1)\n\
		-t t_x t_y			# translation vector (0,0)\n\
		-n N1 N2			# size of output image\n\
		-m mode				# mapping mode (1)\n\
		-c				# keep image centre fixed\n\
		-C				# keep center of (N1/2, N2/2) pixel fixed\n\
		-i				# do inverse transformation\n\
		-f fitsfile			# source for target image\n\
\n\
DESCRIPTION\n\
	\"transformimage\" applies a general linear transformation to a source\n\
	image fs(x) to make a target image f(r) = fs(x(r))\n\
	where the mapping is x_i(r) = psi_ij r_j + t_i.\n\
	By default output image = input image size.\n\
	Use -m option to specify mode, where these are (in order of expense)\n\
		mode = 0:	# nearest pixel\n\
		mode = 1:	# linear interpolation\n\
		mode = 2:	# sum over triangles\n\
	With -c option we calculate t_x, t_y so that the centre\n\
	pixel is mapped to centre of output pixel.\n\
	By default, target image is initialised to zero\n\
	but use -f option to read in an image on which we paint\n\
	the mapped pixels.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../utils/error.h"
#include "../imlib/fits.h"
#include "../imlib/map.h"

int    deflection(float ri, float rj, float *di, float *dj);
float	psi[2][2], t[2];

#define TINY	1.e-50

#define PRESERVE_ORIGIN		0
#define PRESERVE_IMAGE_CENTER	1
#define PRESERVE_PIXEL_CENTER	2

main(int argc, char *argv[])	
{
	int		arg = 1;
	int		N1, N2, M1, M2, mapmode, pixtype, centering, inverse;
	float		**fsource, **ftarget;
	float		e0, e1, a, b, c, d, t0, t1, det;
	FILE		*targetf = NULL;
	fitsheader	*fitsin, *fitsout;
	

	/* defaults */
	psi[0][0] = psi[1][1] = 1.0;
	psi[1][0] = psi[0][1] = 0.0;
	t[0] = t[1] = 0.0;
	M1 = M2 = 0;
	mapmode = FAST_MAP_MODE;
	centering = PRESERVE_ORIGIN;
	inverse = 0;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch (argv[arg++][1]) {
			case 'p':
				if (1 != sscanf(argv[arg++], "%f", &(psi[0][0])) ||
					1 != sscanf(argv[arg++], "%f", &(psi[0][1])) ||
					1 != sscanf(argv[arg++], "%f", &(psi[1][0])) ||
					1 != sscanf(argv[arg++], "%f", &(psi[1][1])))
						error_exit(usage);
				break;
			case 't':
				if (1 != sscanf(argv[arg++], "%f", &(t[0])) ||
					1 != sscanf(argv[arg++], "%f", &(t[1])))
						error_exit(usage);
				break;
			case 'n':
				if (1 != sscanf(argv[arg++], "%d", &M1) ||
					1 != sscanf(argv[arg++], "%d", &M2))
						error_exit(usage);
				break;
			case 'm':
				if (1 != sscanf(argv[arg++], "%d", &mapmode))
						error_exit(usage);
				break;
			case 'c':
				centering = PRESERVE_IMAGE_CENTER;
				break;
			case 'C':
				centering = PRESERVE_PIXEL_CENTER;
				break;
			case 'i':
				inverse = 1;
				break;
			case 'f':
				targetf = fopen(argv[arg++], "r");
				if (!targetf) {
					error_exit("transformimage: failed to open target source\n");
				}
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	if (inverse) {
		a = psi[0][0];
		b = psi[0][1];
		c = psi[1][0];
		d = psi[1][1];
		det = a * d - c * b;
		if (fabs(det) < TINY)
			error_exit("transformimage: zero determinant transformation!\n");
		t0 = t[0];
		t1 = t[1];
		psi[0][0] = d / det;
		psi[0][1] = -b / det;
		psi[1][0] = -c / det;
		psi[1][1] = a / det;
		t[0] = - (psi[0][0] * t0 + psi[0][1] * t1);
		t[1] = - (psi[1][0] * t0 + psi[1][1] * t1);
	}

	/* read the source file */
	read2Dfloatimage(&fsource, &N1, &N2, &fitsin, stdin);
	if (!M1 || !M2) {
		M1 = N1;
		M2 = N2;
	}

	/* create or read the target image */
	if (targetf) {
		read2Dfloatimage(&ftarget, &M1, &M2, &fitsout, targetf);
	} else {
		fitsout = copyfitsheader(fitsin);
		fitsout->n[0] = M1;
		fitsout->n[1] = M2;
		allocFloatArray(&ftarget, M1, M2);
	}
	
	/* subtract delta_ij */
	psi[0][0] -= 1.0;
	psi[1][1] -= 1.0;

	if (centering) {
		t[0] = t[1] = 0;
		switch (centering) {
			case 1:
				deflection(0.5 * N2, 0.5 * N1, &(t[1]), &(t[0]));
				break;
			case 2:
				deflection(0.5 * (N2 + 1), 0.5 * (N1 + 1), &(t[1]), &(t[0]));
				break;
		}
		t[0] *= -1;
		t[1] *= -1;
	}

	/* apply the mapping */
	switch(mapmode) {
		case ULTRAFAST_MAP_MODE:
			ultrafastmap(ftarget, M1, M2, fsource, N1, N2, deflection);
			break;
		case FAST_MAP_MODE:
			fastmap(ftarget, M1, M2, fsource, N1, N2, deflection);
			break;
		case TRIANGLE_MAP_MODE:
			map(ftarget, M1, M2, fsource, N1, N2, deflection);
			break;
		default:
			error_exit("transformimage: bad mapping mode\n");
			break;
	}
	
	
	/* write the target image */
	add_comment(argc, argv, fitsout);
	write2Dfloatimage(ftarget, fitsout);
	
	/* all done */
	exit(0);
}



int	deflection(float y, float x, float *dy, float *dx)
{
	*dx = psi[0][0] * x + psi[0][1] * y + t[0];
	*dy = psi[1][0] * x + psi[1][1] * y + t[1];
	return(1);
}






