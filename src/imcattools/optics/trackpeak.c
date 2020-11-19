#define usage "\n\n\n\
NAME\n\
	trackpeak --- read stream of PSF images and track cell averaged peak\n\
SYNOPSIS\n\
	trackpeak [options...]\n\
\n\
DESCRIPTION\n\
	'trackpeak' reads a 3D FITS image consisting of a stream\n\
	of PSF images, averages them with a cell of size nx * nx * nt, and\n\
	then stacking them with and without shifting.  The purpose is to\n\
	demonstrate the effect of orthogonal transfer devices on PSF quality\n\
\n\
	Options\n\
		-x nx		# cell size in angle (6)\n\
		-t nt		# cell size on time (3)\n\
		-n nphotons	# number of photons (1000)\n\
		-r readnoise	# read noise in electrons (4.0)\n\
		-T		# only output peak position\n\
		-o opfname	# save accumulated psfs as opfname.fits\n\
		-O opfname	# ditto, but don't output the video stream\n\
		-N nframes	# only read the first nframes frames\n\
		-g guidepsf targetpsf	# streams for the guide and target star psfs\n\
		-d dxg dyg dxt dyt	# shift input guide and target star PSFs by this distance\n\
		-G rg		# smooth input PSFs with Gaussian of scale rg\n\
		-f		# do fractional pixel shifting\n\
		-s rs		# smooth the pixellated image before peak finding\n\
		-c		# do centroid guiding\n\
\n\
	In the steady state, the program loops over input guide star image\n\
	planes, shifts each one as it comes in according to shift vectors provided\n\
	with flag -d (zero by default), and installs it in a rotating buffer of\n\
	depth nt.  With -G option it also smooths the image with a Gaussian kernel.\n\
	It computes the average of these planes and installs that in\n\
	another rotating buffer.  Then, from the oldest average, it makes a\n\
	pixellated image and samples this in a poissonian manner and adds\n\
	Gaussian read noise.  It runs a simple peak finder on this object (it\n\
	locates the hottest pixel and then uses a simple shift based on the 1st\n\
	and 2nd derivatives around the hottest pixel.  It accumulates two versions\n\
	of the most recent target images; one unshifted and one shifted according to\n\
	the location of the peak.\n\
\n\
	With -c option we guide on the centroid.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "utils/error.h"
#include "utils/args.h"
#include "imlib/fits.h"
#include "imlib/filters.h"
#include "findpeak.h"
#include "utils/iostream.h"
#include "utils/gaussdev.h"

float 	poidev(float xm, long *idum);
int	getcentroid(float **f, int N1, int N2, double *xc, double *yc);
void    shiftimage(float **fdst, float **fsrc, int Nx, int Ny, int dx, int dy);

int		main(int argc, char *argv[])	
{
	int		x, y, t, nx, nt, N1, N2, p, np, xs, ys, M1, M2, X, Y, ix, iy, ixpk, iypk;
	int		dxg, dyg, dxt, dyt, fractionalshifting, smoothpixels, guideoncentroid;
	int		nphotons, outputposition, outputvideo, useguidestar, gausssmooth;
	float		**fin, ***fb, **ftarget, ***favg, **ftmp, **fshift, **fout, fmax, **facc1, **facc2, **fpix, **fp;
	float		xpk, ypk, fpk, readnoise, sigma;
	double		xc, yc, rg, rs;
	fitsheader	*fitsi, *fitso, *fitsg;
	char		*flag, *psfopfname;
	long		idum;
	FILE		*posopf, *psfopf, *targetf, *guidef;
	iostream	*targetpsfstream, *guidepsfstream;

	/* defaults */
	nx = 6;
	nt = 3;
	nphotons = 1000;
	readnoise = 4.0;
	psfopfname = NULL;
	outputposition = 0;
	np = 0;
	outputvideo = 1;
	useguidestar = 0;
	targetf = stdin;
	gausssmooth = 0;
	fractionalshifting = 0;
	dxt = dyt = dxg = dyg = 0;
	smoothpixels = 0;
	guideoncentroid = 0;

	/* parse args */
	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'x':
				nx = getargi();
				break;
			case 't':
				nt = getargi();
				break;
			case 'n':
				nphotons = getargf();
				break;
			case 'r':
				readnoise = getargf();
				break;
			case 'T':
				outputposition = 1;
				break;
			case 'O':
				outputvideo = 0;
			case 'o':
				psfopfname = getargs();
				psfopf = fopen(psfopfname, "w");
				if (!psfopf) {
					error_exit("trackpeak: cannot open accumulated psf image file\n");
				}
				break;
			case 'N':
				np = getargi();
				break;
			case 'g':
				useguidestar = 1;
				guidepsfstream = openiostream(getargs(), "r");
				guidef = guidepsfstream->f;
				targetpsfstream = openiostream(getargs(), "r");
				targetf = targetpsfstream->f;
				break;
			case 'd':
				dxg = getargi();
				dyg = getargi();
                                dxg = getargi();
                                dyg = getargi();
				break;
			case 'G':
				gausssmooth = 1;
				rg = getargd();
				break;
			case 's':
				smoothpixels = 1;
				rs = getargd();
				break;
			case 'f':
				fractionalshifting = 1;
				break;
			case 'c':
				guideoncentroid = 1;
				break;
			default:
				error_exit(usage);
		}
	}

	/* read the fits header */
	fitsi = readfitsheader(targetf);
	if (fitsi->ndim != 3) {
		error_exit("trackpeak: source image must be 3-dimensional\n");
	}
	N1 = fitsi->n[0];
	N2 = fitsi->n[1];
	if (np) {
		np = (np < fitsi->n[2] ? np : fitsi->n[2]);
	} else {
		np = fitsi->n[2];
	}
	
	if (useguidestar) {
		fitsg = readfitsheader(guidef);
		if (fitsg->ndim != 3) {
			error_exit("trackpeak : guide star image must be 3D\n");
		}
		if (fitsg->n[0] != N1 || fitsg->n[1] != N2) {
			error_exit("trackpeak : guide star image size must match target\n");
		}
	}

	/* size of the pixellated image */
	M1 = N1 / nx;
	M2 = N2 / nx;
	
	/* allocate the buffers nt image planes and averages */	
	fb = (float ***) calloc(nt, sizeof(float **));
	favg = (float ***) calloc(nt, sizeof(float **));
	for (t = 0; t < nt; t++) {
		allocFloatArray(&(fb[t]), N1, N2);
		allocFloatArray(&(favg[t]), N1, N2);
	}

	/* output the fits or cat header as appropriate */
	if (outputposition) {
		posopf = popen("lc -C -n p -N '1 2 xobs' -N '1 2 xinst' -N '1 2 xcent'", "w");
		if (!posopf) {
			error_exit("trackpeak: failed to open output pipe\n");
		}
	} else {
		fitso = copyfitsheader(fitsi);
		fitso->n[0] = 2 * fitsi->n[0];
		fitso->n[1] = 2 * fitsi->n[1];
		fitso->n[2] = np;
		add_comment(argc, argv, fitso);
		if (outputvideo) {
			writefitsheader(fitso);
		}
	}

	/* allocate the cell-averaged image, shifted image */
	allocFloatArray(&fp, N1, N2);
	allocFloatArray(&fshift, N1, N2);
	allocFloatArray(&facc1, N1, N2);
	allocFloatArray(&facc2, N1, N2);
	allocFloatArray(&fout, 2 * N1, 2 * N2);
	allocFloatArray(&ftarget, N1, N2);
	allocFloatArray(&fin, N1, N2);

	/* allocate the pixellated image */
	allocFloatArray(&fpix, M1, M2);

	/* big loop */
	for (p = 0; p < np; p++) {
		/* cycle the buffers */
		ftmp = fb[nt - 1];
		for (t = nt - 1; t > 0; t--) {
			fb[t] = fb[t - 1];
		}
		fb[0] = ftmp;
                ftmp = favg[nt - 1];
                for (t = nt - 1; t > 0; t--) {
                        favg[t] = favg[t - 1];
                }
                favg[0] = ftmp;
		/* read a source image */
		readfitsplane((void *) fin, (useguidestar ? fitsg : fitsi));
		shiftimage(fb[0], fin, N1, N2, dxg, dyg);
		if (gausssmooth) {
			gaussfilter(fb[0], N1, N2, fb[0], rg, rg, 0, 0.0);
		}

		/* take the time average */
		for (y = 0; y < N2; y++) {
			for (x = 0; x < N1; x++) {
				favg[0][y][x] = 0.0;
				for (t = 0; t < nt; t++) {
					favg[0][y][x] += fb[t][y][x];
				}
				favg[0][y][x] /= nt;
			}
		}

		/* generate the pixellated image */
		for (Y = 0; Y < M2; Y++) {
			for (X = 0; X < M1; X++) {
				fpix[Y][X] = 0.0;
				for (y = 0; y < nx; y++) {
					for (x = 0; x < nx; x++) {
						fpix[Y][X] += favg[(p < nt ? p : nt - 1)][nx * Y + y][nx * X + x];
					}
				}
				fpix[Y][X] = poidev(nphotons * fpix[Y][X], &idum) + readnoise * gaussdev();
			}
		}

		/* get peak location */
		if (smoothpixels) {
			gaussfilter(fpix, M1, M2, fpix, rs, rs, 0, 0.0);
		}
		if (guideoncentroid) {
			getcentroid(fpix, M1, M2, &xc, &yc);
			xpk = (float) xc;
			ypk = (float) yc;
		} else {
			findpeak(fpix, M1, M2, &ixpk, &iypk, &xpk, &ypk, &fpk);
		}
		/* scale to source pixel size */
		if (fractionalshifting) {
			xpk *= nx;
			ypk *= nx;
		} else {
			xpk = nx * (0.5 + floor(xpk));
			ypk = nx * (0.5 + floor(ypk));
		}

		/* create the fp image for display */
                for (Y = 0; Y < M2; Y++) {
                        for (X = 0; X < M1; X++) {
                                for (y = 0; y < nx; y++) {
                                        for (x = 0; x < nx; x++) {
                                                fp[nx * Y + y][nx * X + x] = fpix[Y][X] / (nphotons * nx * nx);
                                        }
                                }
                        }
                }

		/* this stuff looks a bit broken */
		if (outputposition) {
			fprintf(posopf, "%5d %14.8lg %14.8lg", p, (double) (xpk), (double ) (ypk));
			xc = xpk;
			yc = ypk;
			if (p > nt / 2) {
				findpeak(fb[nt / 2], N1, N2, &ixpk, &iypk, &xpk, &ypk, &fpk);
				getcentroid(fb[nt / 2], N1, N2, &xc, &yc);
			}
			fprintf(posopf, " %14.8lg %14.8lg %14.8lg %14.8lg\n", (double) (xpk), (double ) (ypk), xc, yc);
			continue;
		}
		ixpk = (int) floor(xpk);
		iypk = (int) floor(ypk);
		/* shift the last image and accumulate */
		if (useguidestar) {
			readfitsplane((void *) fin, fitsi);
		}
		shiftimage(ftarget, fin, N1, N2, dxt, dyt);
                if (gausssmooth) {
                       	gaussfilter(ftarget, N1, N2, ftarget, rg, rg, 0, 0.0);
                }
		for (y = 0; y < N2; y++) {
			for (x = 0; x < N1; x++) {
				xs = (x + ixpk + N1 / 2) % N1;
				ys = (y + iypk + N2 / 2) % N2;
				fshift[y][x] = ftarget[ys][xs];
				facc1[y][x] += ftarget[y][x];
				facc2[y][x] += fshift[y][x];
			}
		}
		if (outputvideo) {
			/* zero out the hot pixel */
			if (ixpk >= 0 && ixpk < N1 && iypk >= 0 && iypk < N2) {
				fp[iypk][ixpk] = 0.0;
			}
			/* construct the output image */
			for (y = 0; y < N2; y++) {
				for (x = 0; x < N1; x++) {
					/* the latest input psf */
					fout[y][x] = fb[0][y][x];
					/* the pixellated image */
					fout[y + N2][x] = fp[y][x];
					/* the unshifted accumulant */
					fout[y][x + N1] = facc1[y][x] / (1 + p);
					/* and the shifted one */
					fout[y + N2][x + N1] = facc2[y][x] / (1 + p);
				}
			}
			if (!outputposition) {
				writefitsplane((void *) fout, fitso);
			}
		}			
	}
	if (psfopfname) {
		/* output the accumulated images */
		fitso->n[0] = N1;
		fitso->n[1] = N2;
		fitso->n[2] = 2;
		fitso->ndim = 3;
		fitso->opstream = psfopf;
		writefitsheader(fitso);
		writefitsplane((void *) facc1, fitso);
		writefitsplane((void *) facc2, fitso);
		writefitstail(fitso);
	}
	exit(0);
}

float ran1(long *idum)
{
	return ((float) drand48());
}


int	getcentroid(float **f, int N1, int N2, double *xc, double *yc)
{
	int	x, y;
	double	fsum, fxsum, fysum;

	fsum = fxsum = fysum = 0.0;
	for (y = 0; y < N2; y++) {
		for (x = 0; x < N1; x++) {
			fsum += f[y][x];
			fxsum += x * f[y][x];
			fysum += y * f[y][x];
		}
	}
	*xc = fxsum / fsum;
	*yc = fysum / fsum;
}

void	shiftimage(float **fdst, float **fsrc, int Nx, int Ny, int dx, int dy)
{
	int x, y, xsrc, ysrc;

	for (y = 0; y < Ny; y++) {
		ysrc = y - dy;
		for (x = 0; x < Nx; x++) {
			xsrc = x - dx;
			if (xsrc >= 0 && xsrc < Nx && ysrc >= 0 && ysrc < Ny) {
				fdst[y][x] = fsrc[ysrc][xsrc];
			} else {
				fdst[y][x] = 0;
			}
		}
	}
}

