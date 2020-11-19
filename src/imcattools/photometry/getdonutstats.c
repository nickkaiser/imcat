#define usage "\nNAME\n\
	getdonutstats - compute shape statistics for out of focus images\n\
\n\
SYNOPSIS\n\
	getdonutstats [-s] [-S] [-m] [-f maskfactor] [-M] [-p nplanes] [-g] [-n nbins (16)] [-N nphi (64)] [-l fsum_min (1e-3)] mmax\n\
\n\
DESCRIPTION\n\
	getdonutstats reads a FITS file containing one or more postage stamp\n\
	images from stdin and computes shape statistics which it sends to stdout\n\
	as a lc-formal calalog.\n\
\n\
ALGORITHM\n\
	For each image segment, getdonutstats first generates an image ff[][] = f[][] * f[][]\n\
	and then calculates the mean of f[][] and of ff[][] for the 'corner' pixels (those\n\
	lying at distance more then n/2 from image center.  It then subtracts these means\n\
	and sets the corner region pixels to zero.  It then computes the flux-weighted centroid\n\
	and computes f0 = sum ff / sumf; r0 = sum fr / sum f; and w0 = (sum f)^2 / 2 pi r0 sum ff\n\
	where r is pixel position relative to the centroid.  It then computes analogous quantities\n\
	for a set of bins in azimuthal angle (16 by default; use -n option to change this) and\n\
	computes the angular Fourier transforms of these for m = 1 through mmax.  The result is\n\
	fm[3][mmax][2] containing, for each of the 3 statistics (flux, radius and width) the real\n\
	and imaginary parts of the transforms.  Radius and width output as is, but the flux is\n\
	multiplied by r0 / f0 to give r0 times the fractional brighness modulation in order\n\
	to obtain a statistic that is independent of defocus distance.\n\
\n\
	With -s option, we read getdonutstats output and synthesise the Fourier\n\
	statistics with nphi steps in azimuth.  Output is a catalog containing 2D image\n\
	index x[2], offset of centroid dx[2], f0, r0, w0 and f[3] containing f, r, w\n\
	and we also output a unit position vector xhat.\n\
\n\
	Option -S is similar to -s, but we just output positions x[2] suitable for\n\
	plotting to show r(phi) and r(phi) +- w(phi)/2.\n\
\n\
	With -m option we output the input image with the regions outside a simple\n\
	estimate of the boundary of the donut set to zero.  Use a maskfactor greater than\n\
	unity to expand the non-masked area.\n\
\n\
	With -M option we output a visualisation of the model as a FITS file.  Similarly, with\n\
	 -p option, except we read the shape statistics catalog (that should contain nplanes\n\
	entries) from stdin.\n\
\n\
	Use -u option to output the man-page.\n\
\n\
SEE ALSO\n\
	sliceimage\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../utils/stats_stuff.h"
#include "../../utils/error.h"
#include "../../utils/arrays.h"
#include "../../imlib/fits.h"
#include "../../utils/args.h"

#define ANALYSIS_MODE	0
#define SYNTHESIS_MODE1	1
#define SYNTHESIS_MODE2 2
#define SYNTHESIS_MODE3 3
#define IPBUFFSIZE	7

int     multiply(double **fm, int mmax, double mulfac);
int     synthesise(double **fm, int mmax, double *fsyn, int nphi);
int     maskimage(float **f, float **r, float **phi, int nx, int ny,  double r0, double w0, double **fsyn, int nphi, double maskfactor);
int	makemodel(float **f, float **r, float **phi, int nx, int ny,  double f0, double r0, double w0, double **fsyn, int nphi, int nstats);

main (int argc, char *argv[]) {
	char		*flag, lccom[1024], tmpstr[1024];
	fitsheader	*fits;
	int		nstats, i, j, nx, ny, iX, iY, ix, iy, m, mmax, opmode, nbins, nphi, iphi, nplanes, dim[3];
	short		**jphi;
	float		**f, **ff, **r, **phi;
	double		nsum, fsum, fxsum, fysum, frsum, ffsum, *fsumphi, *ffsumphi, *frsumphi;
	double		x0, y0, f0, r0, w0, rr, co, si;
	double		***fm, xhat[2], x[2], **fsyn, **fphi, dphi, *ipbuff;
	FILE		*ipstream;
	int		writemaskedimages, makemodelimages;
	double		maskfactor, fsummin;

	/* defaults */
	nbins		= 16;
	nphi		= 64;
	opmode 		= ANALYSIS_MODE;
	writemaskedimages = 0;
	makemodelimages	= 0;
	maskfactor 	= 1.0;
	nstats 		= 3;
	fsummin		= 1.e-3;

        /* parse args */
        argsinit(argc, argv, usage);
        while(nextargtype() == FLAG_ARG) {
		flag = getflag();
                switch (flag[0]) {
			case 's':
				opmode = SYNTHESIS_MODE1;
				break;
			case 'S':
				opmode = SYNTHESIS_MODE2;
				break;
			case 'm':
				writemaskedimages = 1;
				break;
			case 'f':
				maskfactor = getargd();
				break;
			case 'M':
				makemodelimages = 1;
				break;
			case 'n':
				nbins = getargi();
				break;
			case 'N':
				nphi = getargi();
				break;
			case 'p':
				nplanes = getargi();
				opmode = SYNTHESIS_MODE3;
				break;
                        case 'u':
                        default:
                                error_exit(usage);
                }
        }
	if (makemodelimages && writemaskedimages) {
		error_exit("getdonutstats : error : you can't use -m and -M options together\n");
	}
	mmax = getargi();

	dphi = 2 * M_PI / nbins;

        /* allocate statistics */
	fm = (double ***) calloc(nstats, sizeof(double **));
	fphi = (double **) calloc(nstats, sizeof(double *));
	fsyn = (double **) calloc(nstats, sizeof(double *));
        for (i = 0; i < nstats; i++) {
                fm[i] = (double **) calloc(mmax, sizeof(double *));
                for (m = 0; m < mmax; m++) {
                        fm[i][m] = (double *) calloc(2, sizeof(double));
                }
		/* and space for azimuthal stats and synthesised values */
		fphi[i] = (double *) calloc(nbins, sizeof(double));
		fsyn[i] = (double *) calloc(nphi, sizeof(double));
        }
	/* and sums */
	fsumphi = (double *) calloc(nbins, sizeof(double));
	ffsumphi = (double *) calloc(nbins, sizeof(double));
	frsumphi = (double *) calloc(nbins, sizeof(double));

    if (opmode == ANALYSIS_MODE) {
	/* read the header and get image dimensions */
	fits = readfitsheader(stdin);
	if (fits->ndim < 2) {
		error_exit("getdonutstats : error : input image must be at least 2-dimensional");
	}
	nx = fits->n[0];
	ny = fits->n[1];

	/* allocate arrays */
	allocFloatArray(&f, nx, ny);
	allocFloatArray(&ff, nx, ny);
	allocFloatArray(&r, nx, ny);
	allocFloatArray(&phi, nx, ny);

	if (writemaskedimages || makemodelimages) {
		writefitsheader(fits);
	} else {
		/* create the cat header */
		sprintf(lccom, "lc -C -H 'imsize = %d' -N '1 2 x' -N '1 2 dx' -n f0 -n r0 -n w0 -N '3 %d %d 2 fm' < /dev/null", nx, nstats, mmax);
		system(lccom);
	}

	/* loop over image planes */
	for (iY = 0; iY < (fits->n[3] ? fits->n[3] : 1); iY++) {
	    for (iX = 0; iX < (fits->n[2] ? fits->n[2] : 1); iX++) {
		/* read the image */
		readfitsplane((void **) f, fits);
		/* and compute ff[][] = f[][]**2 */
                for (iy = 0; iy < ny; iy++) {
                        for (ix = 0; ix < nx; ix++) {
				ff[iy][ix] = f[iy][ix] * f[iy][ix];
			}
		}
		/* zero the azimuthal sums */
		for (iphi = 0; iphi < nbins; iphi++) {
			fsumphi[iphi] = ffsumphi[iphi] = frsumphi[iphi] = 0.0;
		}
		/* accumulate pixel count, f and f * f in the outer corners (r > nx / 2)  */
		ffsum = fsum = nsum = 0.0;
                for (iy = 0; iy < ny; iy++) {
                        for (ix = 0; ix < nx; ix++) {
				r[iy][ix] = sqrt((iy + 0.5 - ny / 2) * (iy + 0.5 - ny / 2) + (ix + 0.5 - nx / 2) * (ix + 0.5 - nx / 2));
				if (r[iy][ix] > nx / 2) {
					ffsum += ff[iy][ix];
					fsum += f[iy][ix];
					nsum += 1.0;
				}
			}
		}
		/* subtract mean of f and ff computed from outer corners from f[][] and ff[][] images and set outer corners to zero */
                for (iy = 0; iy < ny; iy++) {
                        for (ix = 0; ix < nx; ix++) {
				if (r[iy][ix] < nx / 2) {
					f[iy][ix] -= fsum / nsum;
					ff[iy][ix] -= ffsum / nsum;
				} else {
					f[iy][ix] = ff[iy][ix] = 0.0;
				}
			}
		} 
		/* compute sum f, f * x and f * y */
		fsum = fxsum = fysum = ffsum = frsum = 0.0;
		for (iy = 0; iy < ny; iy++) {
			for (ix = 0; ix < nx; ix++) {
				fsum += f[iy][ix];
				ffsum += ff[iy][ix];
				fxsum += ix * f[iy][ix];
				fysum += iy * f[iy][ix];
			}
		}
		if (fsum > fsummin) {
			/* compute the centroids */
			x0 = fxsum / fsum;
			y0 = fysum / fsum;
			/* accumulate azimuthal sums */
                	for (iy = 0; iy < ny; iy++) {
                	        for (ix = 0; ix < nx; ix++) {
					r[iy][ix] = sqrt((ix - x0) * (ix - x0) + (iy - y0) * (iy - y0));
					frsum += r[iy][ix] * f[iy][ix];
					phi[iy][ix] = atan2(iy - y0, ix - x0);
					if (phi[iy][ix] < 0.0) {
						phi[iy][ix] += 2 * M_PI;
					}
					iphi = (int) floor(phi[iy][ix] / dphi);
					if (iphi == nbins) {
						iphi--;
					}
					fsumphi[iphi] += f[iy][ix];
					ffsumphi[iphi] += ff[iy][ix];
					frsumphi[iphi] += r[iy][ix] * f[iy][ix];
				}
			}
			/* compute mean radius and width */
			r0 = frsum / fsum;
			w0 = fsum * fsum / (2 * M_PI * r0 * ffsum);
			f0 = ffsum / fsum;
			/* compute f(phi), r(phi), w(phi) */
			for (iphi = 0; iphi < nbins; iphi++) {
				fphi[0][iphi] = (((fsumphi[iphi] > 0.0) && (f0 > 0.0))? r0 * ffsumphi[iphi] / (f0 * fsumphi[iphi]) : 0.0);
				fphi[1][iphi] = ((fsumphi[iphi] > 0.0) ? frsumphi[iphi] / fsumphi[iphi] : 0.0);
				fphi[2][iphi] = ((ffsumphi[iphi] > 0.0) ? fsumphi[iphi] * fsumphi[iphi] / ffsumphi[iphi] : 0.0);
				fphi[2][iphi] = ((fphi[1][iphi] > 0.0) ? fphi[2][iphi] / (fphi[1][iphi] * dphi) : 0.0);
			}
			/* initialise Fourier stats */
			for (i = 0; i < nstats; i++) {
				multiply(fm[i], mmax, 0.0);
			}
			/* accumulate */
			for (m = 0; m < mmax; m++) {
				for (iphi = 0; iphi < nbins; iphi++) {
					co = cos((m + 1) * (iphi + 0.5) * dphi);
					si = sin((m + 1) * (iphi + 0.5) * dphi);
					for (i = 0; i < nstats; i++) {
						fm[i][m][0] += fphi[i][iphi] * co;
						fm[i][m][1] += fphi[i][iphi] * si;
					}
				}
			}
			/* normalise the transforms  - with factor 2 because we only had +ve freq - and synthesise */
			for (i = 0; i < nstats; i++) {
				multiply(fm[i], mmax, 2.0 / nbins);
				synthesise(fm[i], mmax, fsyn[i], nphi);
			}
			maskimage(f, r, phi, nx, ny,  r0, w0, fsyn, nphi, maskfactor);
			if (writemaskedimages || makemodelimages) {
				if (makemodelimages) {
					makemodel(f, r, phi, nx, ny,  f0, r0, w0, fsyn, nphi, nstats);
				}
				writefitsplane((void **) f, fits);
			} else {
				/* output */
				printf("%d %d %10.5le %10.5le %10.5le %10.5le %10.5le", iX, iY, x0, y0, f0, r0, w0);
				for (i = 0; i < nstats; i++) {
					for (m = 0; m < mmax; m++) {
                        	        	for (j = 0; j < 2; j++) {
							printf(" %10.5e", fm[i][m][j]);
						}
					}
				}
				printf("\n");
			}
		} else {
			if (writemaskedimages || makemodelimages) {
				writefitsplane((void **) f, fits);
			}
		}
	    }
	}
	} else {
		/* SYNTHESIS_MODE */
		dphi = 2 * M_PI / nphi;
                /* open pipe to read the data */
                ipstream = popen("lc -b -o -P imsize x dx f0 r0 w0 fm", "r");
                /* get the image segment size */
                fgets(tmpstr, 1024, ipstream);
                sscanf(tmpstr, "%d", &nx);
		ny = nx;
		/* generate the output header */
		switch (opmode) {
			case SYNTHESIS_MODE1:
	        		sprintf(lccom, "lc -C -b -N '1 2 x' -N '1 2 dx' -n f0 -n r0 -n w0 -N '1 2 xhat' -N '1 %d f' < /dev/null", nstats);
       				system(lccom);
				break;
			case SYNTHESIS_MODE2:
				system("lc -C -b -N '1 2 x' < /dev/null");
				break;
			case SYNTHESIS_MODE3:
				dim[0] = nx;
				dim[1] = ny;
				dim[2] = nplanes;
				fits = newfitsheader(3, dim, FLOAT_PIXTYPE);
				writefitsheader(fits);
			        allocFloatArray(&f, nx, ny);
			        allocFloatArray(&r, nx, ny);
			        allocFloatArray(&phi, nx, ny);
				break;
			default:
				error_exit("getdonutstats : bad opmode");	
		}
		/* allocate space for x, dx, f0, r0 and w0 */
		ipbuff = (double *) calloc(IPBUFFSIZE, sizeof(double));
		while (IPBUFFSIZE == fread(ipbuff, sizeof(double), IPBUFFSIZE, ipstream)) {
			x0 = ipbuff[2];
			y0 = ipbuff[3];
			f0 = ipbuff[4];
			r0 = ipbuff[5];
			w0 = ipbuff[6];
			for (i = 0; i < nstats; i++) {
				for (m = 0; m < mmax; m++) {
					fread(fm[i][m], sizeof(double), 2, ipstream);
				}
				synthesise(fm[i], mmax, fsyn[i], nphi);
			}
			if (opmode == SYNTHESIS_MODE3) {
	                        for (iy = 0; iy < ny; iy++) {
                                	for (ix = 0; ix < nx; ix++) {
                                        	r[iy][ix] = sqrt((ix - x0) * (ix - x0) + (iy - y0) * (iy - y0));
                                                phi[iy][ix] = atan2(iy - y0, ix - x0);
                                        }
                                }
                                makemodel(f, r, phi, nx, ny, f0, r0, w0, fsyn, nphi, nstats);
                                writefitsplane((void **) f, fits);
			} else {
			    for (iphi = 0; iphi < nphi; iphi++) {
				xhat[0] = cos((iphi + 0.5) * dphi);
				xhat[1] = sin((iphi + 0.5) * dphi);
				switch (opmode) {
				    case SYNTHESIS_MODE1:
					fwrite(ipbuff, sizeof(double), IPBUFFSIZE, stdout);
					fwrite(xhat, sizeof(double), 2, stdout);
					for (i = 0; i < nstats; i++) {
						fwrite(&(fsyn[i][iphi]), sizeof(double), 1, stdout);
					}
					break;
				    case SYNTHESIS_MODE2:
					for (i = 0; i < 2; i++) {
						x[i] = nx * ipbuff[i] + ipbuff[2 + i] + xhat[i] * (ipbuff[5] + fsyn[1][iphi]);
					}
					fwrite(x, sizeof(double), 2, stdout);
                                        for (i = 0; i < 2; i++) {
                                                x[i] = nx * ipbuff[i] + ipbuff[2 + i] + xhat[i] * (ipbuff[5] + fsyn[1][iphi] + 0.5 * (ipbuff[6] + fsyn[2][iphi]));
                                        }
                                        fwrite(x, sizeof(double), 2, stdout);
                                        for (i = 0; i < 2; i++) {
                                                x[i] = nx * ipbuff[i] + ipbuff[2 + i] + xhat[i] * (ipbuff[5] + fsyn[1][iphi] - 0.5 * (ipbuff[6] + fsyn[2][iphi]));
                                        }
                                        fwrite(x, sizeof(double), 2, stdout);
					break;
				}
			    }
			}
		}
	}

	if (writemaskedimages || makemodelimages || (opmode == SYNTHESIS_MODE3)) {
		writefitstail(fits);
	}

	exit(0);
}


int     makemodel(float **f, float **r, float **phi, int nx, int ny,  double f0, double r0, double w0, double **fsyn, int nphi, int nstats)
{
        int     ix, iy, iphi;
        double  dphi, rphi, wphi;

        dphi = 2 * M_PI / nphi;

        for (iy = 0; iy < ny; iy++) {
                for (ix = 0; ix < nx; ix++) {
                        if (phi[iy][ix] < 0.0) {
                        	phi[iy][ix] += 2 * M_PI;
                        }
                        iphi = (int) floor(phi[iy][ix] / dphi);
			if (iphi == nphi) {
				iphi--;
			}
                        rphi = r0 + fsyn[1][iphi];
                        wphi = w0 + fsyn[2][iphi];
                        if ((r[iy][ix] < rphi + 0.5 * wphi) && (r[iy][ix] > rphi - 0.5 * wphi)) {
                                f[iy][ix] = ((r0 > 0.0) ? f0 * (1.0 + fsyn[0][iphi] / r0) : f0);
                        } else {
				f[iy][ix] = 0.0;
			}
                }
        }
}


int	maskimage(float **f, float **r, float **phi, int nx, int ny,  double r0, double w0, double **fsyn, int nphi, double maskfactor)
{
	int 	ix, iy, iphi;
	double	dphi, rphi, wphi;

	dphi = 2 * M_PI / nphi;

	for (iy = 0; iy < ny; iy++) {
		for (ix = 0; ix < nx; ix++) {
			iphi = (int) floor(phi[iy][ix] / dphi);
			rphi = r0 + fsyn[1][iphi];
			wphi = maskfactor * (w0 + fsyn[2][iphi]);
			if ((r[iy][ix] > rphi + 0.5 * wphi) || (r[iy][ix] < rphi - 0.5 * wphi)) {
				f[iy][ix] = 0.0;
			}
		}
	}
}

int	multiply(double **fm, int mmax, double mulfac)
{
	int 	m, j;

	for (m = 0; m < mmax; m++) {
        	for (j = 0; j < 2; j++) {
               		fm[m][j] *= mulfac;
		}
	}
}

int	synthesise(double **fm, int mmax, double *fsyn, int nphi)
{
	int 	i, m;
	double	dphi;

	dphi = 2 * M_PI / nphi;
	for (i = 0; i < nphi; i++) {
		fsyn[i] = 0.0;
                for (m = 0; m < mmax; m++) {
	                fsyn[i] += fm[m][0] * cos((i + 0.5) * (m + 1) * dphi);
                        fsyn[i] += fm[m][1] * sin((i + 0.5) * (m + 1) * dphi);
                }
	}
}

