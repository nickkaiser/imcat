#define usage "\n\n\n\
NAME\n\
	profile - compute azimuthal average of a FITS file\n\
\n\
SYNOPSIS\n\
	profile [option....] < fitsin\n\
		-c xc yc	# spatial origin (N1/2 + 0.5, N2/2 + 0.5)\n\
		-n		# does nothing\n\
		-r rmax		# maximum radius (N1 / 2)\n\
		-d dr		# linear steps - size dr (1.0)\n\
		-l r1 r2 nbins	# logarithmic steps \n\
\n\
DESCRIPTION\n\
	\"profile\" computes the azimuthally averaged profile\n\
	of an image. The output is an lc format catalogue containing\n\
		i r f fsum npix fmin fmax\n\
	where i = 0,1,2...., \n\
	r = i * dr,\n\
	npix is the number of pixels for which the distance\n\
	R from the centre of the pixel to the point (xc,yc)\n\
	lies in the range r <= R < r + 1,\n\
	fsum is the summed image values over those pixels\n\
	and f = fsum / npix.\n\
\n\
	By default the spatial origin is taken to be the centre of\n\
	the pixel (N1/2, N2/2) and the maximum radius\n\
	is half the width of the image: rmax = N1 / 2\n\
\n\
	Use the -l option to do loarithmically spaced bins.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <unistd.h>

#include "../../imlib/fits.h"
#include "../../utils/error.h"

#define      MAX(x,y) (((x) > (y)) ? (x) : (y)) 
#define      MIN(x,y) (((x) < (y)) ? (x) : (y)) 




int		main(int argc, char *argv[])	
{
	int		arg = 1, N1, N2, centreoption;
	fitsheader	*fits;
	float		**f, xc, yc, rmax, dr;
	FILE		*opf;
        float   *nsum, *fsum, *fav, x, y, *fmin, *fmax;
        int     bin, i, j, nbins;
        int     *ntot;
	int	logbins;
	float	r, r1, r2;

	/* defaults */
	centreoption = 0;
	rmax = 0.0;
	dr = 1.0;
	logbins = 0;
	
       	while (arg < argc) {
                if (*argv[arg] != '-')
                        error_exit(usage);
                switch (*(argv[arg++]+1)) {
                        case 'c':
				if (1 != sscanf(argv[arg++], "%f", &xc))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &yc))
					error_exit(usage);
				centreoption = 1;
                                break;
                        case 'r':
				if (1 != sscanf(argv[arg++], "%f", &rmax))
					error_exit(usage);
                                break;
                        case 'd':
				if (1 != sscanf(argv[arg++], "%f", &dr))
					error_exit(usage);
                                break;
                    	case 'n':
                                break;
			case 'l':
				logbins = 1;
				if (1 != sscanf(argv[arg++], "%f", &r1))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &r2))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%d", &nbins))
					error_exit(usage);
				break;
                        default:
                                error_exit(usage);
                                break;
                }
        }
        
	opf = popen("lc -C -n i -n r -n f -n fsum -n npix -n fmin -n fmax", "w");
	if (!opf) {
		fprintf(stderr, "profile: failed to open output pipe\n");
		exit(-1);
	}

	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);

	if (!centreoption) {
		xc = N1 / 2 + 0.5;
		yc = N2 / 2 + 0.5;
	}

	if (rmax == 0.0) {
		rmax = 0.5 * N1;
	}

         /* set up the bins */
	if (!logbins) {
	        nbins  = ceil(rmax / dr);
	}
        nsum = (float *) calloc(nbins, sizeof(float));
        fsum = (float *) calloc(nbins, sizeof(float));
        fav = (float *) calloc(nbins, sizeof(float));
        fmin = (float *) calloc(nbins, sizeof(float));
        fmax = (float *) calloc(nbins, sizeof(float));
 	for (i = 0; i < nbins; i++) {
		fmax[i] = -FLT_MAX;
		fmin[i] =  FLT_MAX;
	}

        /* accumulate nsum, fsum. ffsum */
        for (i = 0; i < N2; i++) {
		y = i + 0.5;
                for (j = 0; j < N1; j++) {
			x = j + 0.5;
			r = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
			if (logbins) {
				if (r > 0.0) {
					bin = (int) floor(nbins * log(r / r1) / log(r2 / r1));
				}  
			} else {
	                        bin = (int) floor(r / dr);
			}
                        if (f[i][j] == FLOAT_MAGIC)
                                continue;
                        if (bin < 0 || bin >= nbins)
                                continue;
                        nsum[bin] += 1.0;
                        fsum[bin] += (float) f[i][j];
			if (f[i][j] < fmin[bin]) {
				fmin[bin] = f[i][j];
			}  
			if (f[i][j] > fmax[bin]) {
				fmax[bin] = f[i][j];
			}  
                }
        }
        for (bin = 0; bin < nbins; bin++) {
                if (nsum[bin]) {
                        fav[bin] = fsum[bin] / nsum[bin];
                }
        }

        for (bin = 0; bin < nbins; bin++) {
		if (logbins) {
			r = r1 * exp((bin + 0.5) * log(r2 / r1) / nbins);
		} else {
			r = dr * (bin + 0.5);
		}
		fprintf(opf, "%10d %13.8g %13.8g %13.8g %13.8g %13.8g %13.8g\n", 
			bin, r, fav[bin], fsum[bin], nsum[bin], fmin[bin], fmax[bin]);
	}
	exit(pclose(opf));
}



