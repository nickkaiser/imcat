/*
 * warpimage.c
 */

#define	usage "\n\n\n\
NAME\n\
	warpimage - apply spatial transformation to a FITS image\n\
\n\
SYNOPSIS\n\
	warpimage [options...] distparfile\n\
		-n N1 N2	# size of output image\n\
		-g sf N1 N2 	# generate distortion image\n\
		-m		# initialise image to magic\n\
		-M mode		# mapping mode (1)\n\
		-q qparfile	# recircularise\n\
\n\
DESCRIPTION\n\
	\"warpimage\" reads a fits file from stdin and applies a\n\
	spatial transformation according\n\
	to the parameters in 'distparfile' such that:\n\
		r = x + sum_m a_m f_m(x)\n\
	where mode functions are polynomials in x[0], x[1].\n\
	Use -p and -d options to apply a further linear transformation:\n\
		r_i => phi_ij r_j + d_j\n\
	By default output image = input image size.\n\
	Output image is initialised to zero, unless you use -m option.\n\
\n\
	Use -M option to specify mode, where these are (in order of expense)\n\
		mode = 0:	# nearest pixel\n\
		mode = 1:	# linear interpolation\n\
		mode = 2:	# sum over triangles   \n\
\n\
	With -g option (and s = 0) we generate a N1 by (2 * N2) image\n\
	whose lower and upper halves contain the x, y components\n\
	of the distortion d = r-x. Set 'sf' to be some small\n\
	integer to create a deflection image which has been demagnified\n\
	by a factor 2^sf, but N1, N2 are given in source image pixel\n\
	scale, so '-g 3 1024 1024' will generate a 128 x 256 pixel image.\n\
\n\
	With -q option we read a parameter file for a model of the field\n\
	q = e / P_sm and then make appropriate additional deflections\n\
	to recircularise the psf.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@ifa.hawaii.edu\n\
\n\n\n"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../utils/error.h"
#include "../imlib/map.h"
#include "../imlib/fits.h"
#include "../utils/modefunc.h"

void    	deflection(float ri, float rj, float *di, float *dj);
void    	recircdef(float ri, float rj, float *di, float *dj);
int		nmodes, *l, *m;
double		**a;
int		nqmodes, *lq, *mq, defsign;
double		**qpar;
double		scalefac;

main(int argc, char *argv[])	
{
	int		arg = 1;
	char		line[1024], *qparfile;
	FILE		*distparfile;
	int		mode;
	int		dodisplacement, dophitrans, makedistortionimage;
	int		N1, N2, M1, M2, pixtype, ix, iy;
	float		**fsource, **ftarget, **ft1, **ft2, **dx, **dy, di, dj;
	int		scrunchfac, ndim, nqdim, dim[3];
	char		*vardef[MODEFUNC_MAX_VARS], *xname;
	int		nvar, magicinit, mapmode, recircularisepsf;
	fitsheader	*fits;

	/* defaults */
	dodisplacement = 0;
	dophitrans = 0;
	makedistortionimage = 0;
	scrunchfac = 0;
	scalefac = 1.0;
	magicinit = 0;
	M1 = M2 = 0;
	mapmode = FAST_MAP_MODE;
	recircularisepsf = 0;

	/* parse args */
	if (argc < 2)
		error_exit(usage);
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			break;
		} else {
			switch (argv[arg++][1]) {
				case 'g':
					makedistortionimage = 1;
					if ((argc - arg) < 1)
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%d", &scrunchfac))
						error_exit(usage);
				case 'n':
					if ((argc - arg) < 2)
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%d", &M1) ||
						1 != sscanf(argv[arg++], "%d", &M2))
						error_exit(usage);
					break;
				case 'm':
					magicinit = 1;
					break;
				case 'M':
					if ((argc - arg) < 1)
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%d", &mapmode))
						error_exit(usage);
					break;
				case 'q':
					recircularisepsf = 1;
					qparfile = argv[arg++];
					break;
				case 'u':
				default:
					error_exit(usage);
					break;
			}
		}
	}
	if (recircularisepsf) {
		if (mapmode != 1) {
			error_exit("warpimage: you must use mapmode = 1 with -q option\n");
		}
		if (makedistortionimage) {
			error_exit("warpimage: you can't use -q and -g options together\n");
		}
		get2Dpolymodel(qparfile, &lq, &mq, &nqdim, &qpar, &nqmodes, &nvar, vardef, &xname);
		if (nqdim != 2) {
			error_exit("warpcat: I need a 2D q-parameter file!\n");
		}		
	}

	get2Dpolymodel(argv[arg], &l, &m, &ndim, &a, &nmodes, &nvar, vardef, &xname);
	if (ndim != 2) {
		error_exit("warpcat: I need a 2D parameter file!\n");
	}

	
	if (makedistortionimage) {
		while (scrunchfac-- > 0) {
			M1 /= 2;
			M2 /= 2;
			scalefac *= 2;
		}
		dim[0] = M1;
		dim[1] = M2;
		dim[2] = 2;
		fits = newfitsheader(3, dim, FLOAT_PIXTYPE);
		allocFloatArray(&dx, M1, M2);
		allocFloatArray(&dy, M1, M2);
		for (iy = 0; iy < M2; iy++) {
			for (ix = 0; ix < M1; ix++) {
				deflection((float) (scalefac * iy), (float) (scalefac * ix), &di, &dj);
				dx[iy][ix] = dj;
				dy[iy][ix] = di;
			}
		}
		add_comment(argc, argv, fits);
		writefitsheader(fits);
		writefitsplane((void **) dx, fits);
		writefitsplane((void **) dy, fits);
		writefitstail(fits);
		exit(0);
	}

	read2Dfloatimage(&fsource, &N1, &N2, &fits, stdin);
	/* set the target file size if necessary */
	if (!M1 || !M2) {
		M1 = N1;
		M2 = N2;
	}
	if (recircularisepsf) {
		allocFloatArray(&ft1, M1, M2);
		allocFloatArray(&ft2, M1, M2);
		if (magicinit) {
			for (iy = 0; iy < M2; iy++) {
				for (ix = 0; ix < M1; ix++) {
					ft1[iy][ix] = ft2[iy][ix] = FLOAT_MAGIC;
				}
			}
		}
	} else {
		allocFloatArray(&ftarget, M1, M2);
		if (magicinit) {
			for (iy = 0; iy < M2; iy++) {
				for (ix = 0; ix < M1; ix++) {
					ftarget[iy][ix] = FLOAT_MAGIC;
				}
			}
		}
	}
	/* apply the mapping */
	switch(mapmode) {
		case ULTRAFAST_MAP_MODE:
			ultrafastmap(ftarget, M1, M2, fsource, N1, N2, deflection);
			break;
			case FAST_MAP_MODE:
			if (recircularisepsf) {
				defsign = 1;
				fastmap(ft1, M1, M2, fsource, N1, N2, recircdef);
				defsign = -1;
				fastmap(ft2, M1, M2, fsource, N1, N2, recircdef);
				for (iy = 0; iy < M2; iy++) {
					for (ix = 0; ix < M1; ix++) {
						if ((ft1[iy][ix] == FLOAT_MAGIC) || ft2[iy][ix] == FLOAT_MAGIC) {
							ft1[iy][ix] = FLOAT_MAGIC;
						} else {
							ft1[iy][ix] = 0.5 * (ft1[iy][ix] + ft2[iy][ix]);
						}
					}
				}
				ftarget = ft1;
			} else {
				fastmap(ftarget, M1, M2, fsource, N1, N2, deflection);
			}
			break;
		case TRIANGLE_MAP_MODE:
			map(ftarget, M1, M2, fsource, N1, N2, deflection);
			break;
		default:
			error_exit("transformimage: bad mapping mode\n");
			break;
	}

	/* write the target image */
	add_comment(argc, argv, fits);
	fits->n[0] = M1;
	fits->n[1] = M2;
	write2Dfloatimage(ftarget, fits);
	
	/* all done */
	exit(0);
}


void    	deflection(float ri, float rj, float *di, float *dj)
{
	int	i, mode;
	double	fmode, r[2], x[2];

	x[0] = (double) rj;
	x[1] = (double) ri;
	for (i = 0; i < 2; i++) {
		r[i] = 0.0;
	}
	for (mode = 0; mode < nmodes; mode++) {
		fmode = f(l[mode], m[mode], x);
		for (i = 0; i < 2; i++) {
			r[i] += a[i][mode] * fmode;
		}
	}
	*dj = (float) r[0];
	*di = (float) r[1];
}

void    	recircdef(float ri, float rj, float *di, float *dj)
{
	int	i, mode;
	double	fmode, r[2], x[2], q[2], phi, d[2], dmod;

	/* initialise x[] */
	x[0] = (double) rj;
	x[1] = (double) ri;

	/* compute the regular source coordinate */
	for (i = 0; i < 2; i++) {
		r[i] = x[i];
	}
	for (mode = 0; mode < nmodes; mode++) {
		fmode = f(l[mode], m[mode], x);
		for (i = 0; i < 2; i++) {
			r[i] += a[i][mode] * fmode;
		}
	}
	/* compute the q-vector */
	for (i = 0; i < 2; i++) {
		q[i] = 0.0;
	}
	for (mode = 0; mode < nqmodes; mode++) {
		fmode = f(lq[mode], mq[mode], r);
		for (i = 0; i < 2; i++) {
			q[i] += qpar[i][mode] * fmode;
		}
	}
	/* compute the extra deflection */
	phi = 0.5 * (atan2(q[1], q[0]) + M_PI);
	dmod = pow(q[0] * q[0] + q[1] * q[1], 0.25);
	d[0] = defsign * dmod * cos(phi);
	d[1] = defsign * dmod * sin(phi);

	*dj = (float) (r[0] + d[0] - x[0]);
	*di = (float) (r[1] + d[1] - x[1]);
}

