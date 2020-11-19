/**
 **
 ** calculate simple statistics for an image
 **
 **/

#define	usage "\n\n\n\
NAME\n\
	linesmode --- calculate mode of rows/cols of a fits file\n\
\n\
SYNOPSIS\n\
	linesmode [options...]\n\
		-c			# work in columns mode\n\
		-f fitsf		# fix fits file 'fitsf'\n\
\n\
DESCRIPTION\n\
	By default, \"linesmode\"  reads a N1 * N2 fits file from stdin,\n\
	calculates the median for each row (might be useful to\n\
	zap detected objects with 'makechart' first) and writes the\n\
	result to stdout as a 1-D fits file.\n\
	With -c option, we calculate the median for each column instead.\n\
	With -f option we read the line modes (or column modes with -c) mode\n\
	file from stdin, and subtract mode from the named fits\n\
	file, writing the result to standard out.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"		

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../utils/stats_stuff.h"
#include "../utils/error.h"
#include "../imlib/fits.h"

#define LINE_ORIENT_MODE	0
#define COLS_ORIENT_MODE	1
#define GET_MODE_MODE		0
#define FIX_FITS_MODE		1

#define MIN_SAMPLE_SIZE		10

#define MAGIC	FLOAT_MAGIC	

int main(int argc, char *argv[])	
{
	int		arg = 1, N1, N2, pixtype, opmode, orientmode, x, y, nsample;
	float		**f, *fsample, *fout, fmode, median, uq, lq, sigma, *modef;
	fitsheader	*fits, *modefits;
	char		*fitsfilename;
	FILE		*fitsf;


	/* defaults */
	orientmode = LINE_ORIENT_MODE;
	opmode = GET_MODE_MODE;
	fitsf = stdin;
	
	while (arg < argc) {
		if (*argv[arg] != '-') {
			error_exit(usage);
		}
		switch (*(argv[arg++]+1)) {
			case 'f':
				fitsfilename = argv[arg++];
				opmode = FIX_FITS_MODE;
				if (!(fitsf = fopen(fitsfilename, "r")))
					error_exit("linesmode: unable to open fits file\n");
				break;
			case 'c':
				orientmode = COLS_ORIENT_MODE;
				break;
			default:
				error_exit(usage);
				break;
		}
	}
	
	switch (opmode) {
		case GET_MODE_MODE:
			read2Dfloatimage(&f, &N1, &N2, &fits, fitsf);
			fits->ndim = 1;
			fits->extpixtype = FLOAT_PIXTYPE;
			switch (orientmode) {
				case LINE_ORIENT_MODE:
					fsample = (float *) calloc(N1, sizeof(float));
					fout = (float *) calloc(N2, sizeof(float));
					fits->n[0] = N2;
					for (y = 0; y < N2; y++) {
						nsample = 0;
						for (x = 0; x < N1; x++) {
							if (f[y][x] != MAGIC) {
								fsample[nsample++] = f[y][x];
							}
						}
						if (nsample >= MIN_SAMPLE_SIZE) {
							liststats(fsample, nsample, &median, &lq, &uq, &sigma);
							fout[y] = median;
						} else {
							fout[y] = (float) MAGIC;
						}
					}
 					break;
				case COLS_ORIENT_MODE:
					fsample = (float *) calloc(N2, sizeof(float));
					fout = (float *) calloc(N1, sizeof(float));
					fits->n[0] = N1;
					for (x = 0; x < N1; x++) {
						nsample = 0;
						for (y = 0; y < N2; y++) {
							if (f[y][x] != MAGIC) {
								fsample[nsample++] = f[y][x];
							}
						}
						if (nsample >= MIN_SAMPLE_SIZE) {
							liststats(fsample, nsample, &median, &lq, &uq, &sigma);
							fout[x] = median;
						} else {
							fout[x] = (float) MAGIC;
						}
					}
 					break;
				default:
					error_exit("linesmode: bad orientmode\n");
					break;
			}
			writefitsheader(fits);
			writefitsline(fout, fits);
			break;
		case FIX_FITS_MODE:
			modefits = readfitsheader(stdin);
			if (modefits->ndim != 1) {
				error_exit("linesmode : expecting a 1D FITS image on stdin\n");
			}
			modef = (float *) calloc(fits->n[0], sizeof(float));
			readfitsline(modef, modefits);
			read2Dfloatimage(&f, &N1, &N2, &fits, fitsf);
			add_comment(argc, argv, fits);
			switch (orientmode) {
				case LINE_ORIENT_MODE:
					for (y = 0; y < N2; y++) {
						fscanf(stdin, "%f", &fmode);
						for (x = 0; x < N1; x++) {
							if (f[y][x] != MAGIC)
								f[y][x] -= modef[y];
						}
					}
					break;
				case COLS_ORIENT_MODE:
 					for (x = 0; x < N1; x++) {
						fscanf(stdin, "%f", &fmode);
						for (y = 0; y < N2; y++) {
							if (f[y][x] != MAGIC)
								f[y][x] -= modef[x];
						}
					}
					break;
				default:
					error_exit("linesmode: bad orientmode\n");
					break;
			}
			write2Dfloatimage(f, fits);
			break;
		default:
			error_exit("linesmode: bad opmode!\n");
			break;
	}
	exit(0);

}





