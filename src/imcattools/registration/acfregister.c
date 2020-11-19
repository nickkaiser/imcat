 /*
 * acfregister.c
 *
 * determine rotation, scale and offset for pair of lists of coords
 */


#include <stdio.h>
#include <math.h>
#include "utils/arrays.h"
#include "fftlib/myfft.h"
#include "utils/error.h"
#include "utils/ipbuff.h"
#include "imlib/fits.h"
#include "imlib/filters.h"

#define usage "\n\n\
NAME\n\
	acfregister - determine approximate transformation coefficients\n\
\n\
SYNOPSIS\n\
	acfregister [options...] a.cat b.cat\n\
		-x xname	# name for spatial coords ('x')\n\
		-i imsize	# internal image size (512)\n\
		-p phi1 phi2	# range for phi wrapping (0, PI)\n\
		-a a1 a2	# die if solution has scale factor outside this range\n\
		-d dlnr		# range of (natural) log(r) for wrapping (1.0)\n\
		-v 		# pipe images to 'iis -p ~/dev/saopipe' for display\n\
DESCRIPTION\n\
	'acfregister' reads 2 catalogues and determines\n\
	scale, rotation and translation that maps (x_a,y_a) => (x_b,y_b)\n\
\n\
		x_b = a (x_a cos phi - y_a sin phi) + x0\t\n\
		y_b = a (x_a sin phi + y_a cos phi) + y0\t\n\
\n\
	It works by cross-correlating images of wrapped lnr, phi\n\
	values for pairs of objects and locating peak to determine\n\
	a, phi.  We then rotate and scale x_b y_b, generate images\n\
	of x,y positions and cross correlate these to get x0, y0.\n\
\n\
	Outputs x0, y0, a, phi to stdout.\n\
\n\
	Use the -v option to get some visual feedback.  You must first\n\
	have saoimage running and listening to the FIFO ~/dev/saopipe\n\
	so make this pipe if necessary (using mknod) and then run the\n\
	image display using 'saoimage -idev ~/dev/saopipe'.  Acfregister\n\
	will then pipe the various images it generates internally to the viewer.\n\
\n\
	Acfregister will not work properly if the scale difference is\n\
	very large (|log(a)| > dlnr/2) as the wrapping will cause it to\n\
	return log(a) in range +- dlnr/2.  You can increase dlnr, but at\n\
	the cost of reduced precision in log(a).  It is probably better to\n\
	rescale one of the input cats to get roughly similar coordinate scales.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n"

#define	PI M_PI
#define BIG_NEG		-1.e20
#define BIG_POS		1.e20;
#define TINY    	1.e-20

#define DISPLAY_COMMAND "iis -p ~/dev/saopipe"

float   mexicanfilterfunction(float ki, float kj);
void	display(float **f, int N1, int N2, char *name);
int	findpeak(float **f, int imsize, double *xpeak, double *ypeak, double* peakval);

main(int argc, char *argv[])
{
	char	line[1024], defxname[2] = "x", *xname, lccommand[64];
	double	rr, dx[2], dy[2], Phi[2], dphi, phi1, phi2, lnr, dlnr, a1, a2, temp, peakval[2], junk, Dx, Dy;
	int	arg, i, ix, iy, imsize, nobj[2], obj, obj0, obj1, displayimages, cat, iphi; 
	FILE	*ipf;
	float	**n[2], **theccf, sigma1, sigma2;
	double	phi, phimax, amax, c, s, xrange, xpeak, ypeak;
	double	**x[2], **X[2];			/* x[cat][obj][coord] and (transformed) copy */
	double	xll[2][2], xur[2][2];	/* bounding boxes xll[cat][coord] */
	fft_type	fk[2];

	/* defaults */
	phi1 = 0.0;
	phi2 = PI;
	dlnr = 1.0;
	imsize = 512;
	a1 = 0.1;
	a2 = 10.0;
	sigma1 = 1;
	sigma2 = 2;
	xname = defxname;
	displayimages = 0;

	if (argc < 3) {
		fprintf(stderr, "%s", usage);
		exit(-1);
	}

	/* parse args */
	arg = 1;
	while (arg < argc) {
		if (argv[arg][0] != '-')
			break;
		switch (argv[arg++][1]) {
			case 'x':
				xname = argv[arg++];
				break;
			case 'i':
				sscanf(argv[arg++], "%d", &imsize);
				break;
			case 'p':
				sscanf(argv[arg++], "%lf", &phi1);
				sscanf(argv[arg++], "%lf", &phi2);
				break;
			case 'a':
				sscanf(argv[arg++], "%lf", &a1);
				sscanf(argv[arg++], "%lf", &a2);
				if (a2 < a1) {
					temp = a2;
					a2 = a1;
					a1 = temp;
				}
				break;
			case 'd':
				sscanf(argv[arg++], "%lf", &dlnr);
				break;
			case 'v':
				displayimages = 1;
				break;
			default:
				error_exit(usage);
				break;
		}
	}
	dphi = phi2 - phi1;


	/* open the catalogues and read objects into expandable buffers */
	if (arg > argc - 2)
		error_exit(usage);
	for (cat = 0; cat < 2; cat++) {
		sprintf(lccommand, "lc -b -o %s < %s", xname, argv[arg++]);
		if (!(ipf = popen(lccommand,"r")))
			error_exit("register: failed to open pipe for input\n");
		x[cat] = readdoublebuff(2, ipf, &(nobj[cat]));
		pclose(ipf);
		/* allocate space for transformed and shifted versions of cats */
		X[cat] = (double **) calloc(nobj[cat], sizeof(double *));
		for (obj = 0; obj < nobj[cat]; obj++) {
			X[cat][obj] = (double *) calloc(2, sizeof(double));
		}
	}


	/* allocate image arrays, fft's */
	allocFloatArray(&theccf, imsize, imsize);
	for (i = 0; i < 2; i++) {
		allocFloatArray(&(n[i]), imsize, imsize);
		alloc_fft(&(fk[i]), imsize, imsize);
	}

	/* make images of pairs */
	for (cat = 0; cat < 2; cat++) {
		for (obj0 = 0; obj0 < nobj[cat]; obj0++) {
			for (obj1 = 0; obj1 < nobj[cat]; obj1++) {
				Dx = x[cat][obj1][0] - x[cat][obj0][0];
				Dy = x[cat][obj1][1] - x[cat][obj0][1];
				rr = Dx * Dx + Dy * Dy;
				if (rr == 0.0)
					continue;
				lnr = 0.5 * log(rr);
				phi = atan2(Dy, Dx) - phi1;
				ix = imsize * (phi/dphi - floor(phi/dphi));
				iy = imsize * (lnr/dlnr - floor(lnr/dlnr));
				n[cat][iy][ix] += 1.0;
			}
		}
	}

	/* do the crosscorrelation */
	if (displayimages) display(n[0], imsize, imsize, "pairs n[0]");
	forward_fft(n[0], imsize, imsize, fk[0]);
	if (displayimages) display(n[1], imsize, imsize, "pairs n[1]");
	forward_fft(n[1], imsize, imsize, fk[1]);
	ccf(fk[0], fk[1], imsize, imsize, theccf, imsize / 2, imsize / 2);
	if (displayimages) display(theccf, imsize, imsize, "scale-angle-ccf");
	mexicanfilter(theccf, imsize, imsize, theccf, sigma1, sigma2, 0);
	if (displayimages) display(theccf, imsize, imsize, "smoothed scale-angle-ccf");

	/* find the cross correlation peak */
	findpeak(theccf, imsize, &xpeak, &ypeak, &junk);
	amax = exp(dlnr * (ypeak - imsize / 2) / imsize);
	phimax = phi1 + dphi * (xpeak - imsize / 2) / imsize; 
	if (amax < a1 || amax > a2) {
		error_exit("acfregister: peak outside of range a1-a2\n");
	}
	

	/* for each possible orientation... */
	Phi[0] = phimax;
	Phi[1] = phimax + M_PI;
	for (iphi = 0; iphi < 2; iphi++) {
		/* rotate the 1st list of coords */
		c = cos(Phi[iphi]);
		s = sin(Phi[iphi]);
		for (obj = 0; obj < nobj[0]; obj++) {
			X[0][obj][0] = amax * (c * x[0][obj][0] - s * x[0][obj][1]);
			X[0][obj][1] = amax * (s * x[0][obj][0] + c * x[0][obj][1]);
		}
		/* and make copy of second */	
		for (obj = 0; obj < nobj[1]; obj++) {
			X[1][obj][0] = x[1][obj][0];
			X[1][obj][1] = x[1][obj][1];
		}
		/* for each cat... */
		for (cat = 0; cat < 2; cat++) {
			/* figure out the bounding boxes */
			for (i = 0; i < 2; i++) {
				xll[cat][i] = BIG_POS;
				xur[cat][i] = BIG_NEG;
				for (obj = 0; obj < nobj[cat]; obj++) {
					if (X[cat][obj][i] < xll[cat][i]) {
						xll[cat][i] = X[cat][obj][i];
					}
					if (X[cat][obj][i] > xur[cat][i]) {
						xur[cat][i] = X[cat][obj][i];
					}
				}
			}
			/* and shift the origin to xll */
			for (obj = 0; obj < nobj[cat]; obj++) {
				for (i = 0; i < 2; i++) {
					X[cat][obj][i] -= xll[cat][i];
				}
			}
		}
		/* figure out the biggest x-range */
		xrange = 0;
		for (cat = 0; cat < 2; cat++) {
			for (i = 0; i < 2; i++) {
				if ((xur[cat][i] - xll[cat][i]) > xrange) {
					xrange = xur[cat][i] - xll[cat][i];
				}
			}
		}
	
		/* now generate the two images of x,y coords */
		for (cat = 0; cat < 2; cat++) {
			for (iy = 0; iy < imsize; iy++) {
				for (ix = 0; ix < imsize; ix++) {
					n[cat][iy][ix] = 0.0;
				}
			}
			for (obj = 0; obj < nobj[cat]; obj++) {
				ix = floor(imsize * X[cat][obj][0] / xrange);
				iy = floor(imsize * X[cat][obj][1] / xrange);
				if (ix >= 0 && ix < imsize && iy >= 0 && iy < imsize)
					n[cat][iy][ix] += 1.0;
			}
		}

		/* now we cross correlate and smooth */
		if (displayimages) display(n[0], imsize, imsize, "scaled-rotated n[0]");
		forward_fft(n[0], imsize, imsize, fk[0]);
		if (displayimages) display(n[1], imsize, imsize, "scaled-rotated n[1]");
		forward_fft(n[1], imsize, imsize, fk[1]);
		ccf(fk[0], fk[1], imsize, imsize, theccf, imsize / 2, imsize / 2);
		if (displayimages) display(theccf, imsize, imsize, "x-y-ccf");
		mexicanfilter(theccf, imsize, imsize, theccf, sigma1, sigma2, 0);
		if (displayimages) display(theccf, imsize, imsize, "smoothed x-y-ccf");

		/* and find the peak */
		findpeak(theccf, imsize, &xpeak, &ypeak, &(peakval[iphi]));
		/* compute shift */
		dx[iphi] = xrange * (xpeak - imsize / 2) / imsize + xll[1][0] - xll[0][0]; 
		dy[iphi] = xrange * (ypeak - imsize / 2) / imsize + xll[1][1] - xll[0][1];
	}	

	/* choose the highest peak */
	iphi = (peakval[0] > peakval[1] ? 0 : 1);

	/* and we're all done
	fprintf(stdout, "#       dx         dy          a        phi\n"); */
	fprintf(stdout, "%13.8lg %13.8lg %13.8lg %13.8lg\n", dx[iphi], dy[iphi], amax, Phi[iphi]);
	exit(0);	
}

void	display(float **f, int N1, int N2, char *name)
{
	char	*argv[1];
 	fitsheader	*fits;

	fprintf(stderr, "# acfregister: displaying %s   hit return to continue\n", name);
	fits = new2Dfitsheader(N1, N2, FLOAT_PIXTYPE);
	if (!(fits->opstream = popen(DISPLAY_COMMAND, "w")))
		error_exit("register: failed to open display pipe\n");
	write2Dfloatimage(f, fits);
	pclose(fits->opstream);
	getchar();
}


int	findpeak(float **f, int imsize, double *xpeak, double *ypeak, double *thepeakval)
{
	int	i, j, i0, j0, ixpeak, iypeak;	
	double	Fp, Fpp, peakval = BIG_NEG;
	static 	float **ff = NULL;

	/* allocate space for 3*3 array around peak if necessary */
	if (!ff) {
		ff = (float **) calloc(3, sizeof(float *)) + 1;
		for (i0 = -1; i0 <= 1; i0++) {
			ff[i0] = (float *) calloc(3, sizeof(float)) + 1;
		}
	}

	/* find hottest pixel */
	for (i = 0; i < imsize; i++) {
		for (j = 0; j < imsize; j++) {
			if (f[i][j] > peakval) {
				iypeak = i;
				ixpeak = j;
				peakval = f[i][j];
			} 
		}
	}
	*xpeak = (double) ixpeak;
	*ypeak = (double) iypeak;
	
	/* copy 3 x 3 array */
	for (i0 = -1; i0 <= 1; i0++) {
		i = iypeak + i0;
		if (i < 0) {
			i += imsize;
		}
		if (i >= imsize) {
			i -= imsize;
		}
		for (j0 = -1; j0 <= 1; j0++) {
			j = ixpeak + j0;
			if (j < 0) {
				j += imsize;
			}
			if (j >= imsize) {
				j -= imsize;
			}
			ff[i0][j0] = f[i][j];
		}
	}
	
	/* now we add the simple cludge to get more precise position */
	Fp  = 0.5 * (ff[0][1] - ff[0][-1]);
	Fpp = ff[0][1] - 2 * ff[0][0] + ff[0][-1];
	if (fabs(Fpp) > TINY)
		*xpeak -= Fp / Fpp;
	Fp  = 0.5 * (ff[1][0] - ff[-1][0]);
	Fpp = ff[1][0] - 2 * ff[0][0] + ff[-1][0];
	if (fabs(Fpp) > TINY)
		*ypeak -= Fp / Fpp;
	*thepeakval = peakval;
}