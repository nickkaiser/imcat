#define usage "\n\n\n\
NAME\n\
	album --- paste fits images into a larger image\n\
\n\
SYNOPSIS\n\
	album [options... ] nx ny fits1 fits2 ...\n\
	album -f offsetsfile [options... ] N1 N2 fits1 fits2 ...\n\
		-d dx dy	# grid spacing\n\
		-b		# draw a line around each chip\n\
		-f file		# file for offsets\n\
		-M		# initialise to MAGIC.\n\
\n\
DESCRIPTION\n\
	\"album\" combines a set of images into an album.\n\
	With -f option the layout of the image is determined\n\
	from the'offsetsfile', amd the size of the output\n\
	image is specified after the options.\n\
	The format of this file should be a single comment line\n\
	followed by x,y pairs to define location of bottom left\n\
	corners of the images to be pasted.\n\
	Without the -f option, the grid spacing will be\n\
	equal to the dimensions of the first image (unless\n\
	you override this with -d option) and the images will\n\
	be placed on a nx by ny grid.\n\
\n\
	The output image is initialised to zero, images are painted\n\
	on sequentially, erasing any previously painted values, except\n\
	that MAGIC values are not painted.  Us the -M option to\n\
	intialise to MAGIC instead.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/arrays.h"


int		main(int argc, char *argv[])	
{
	int	arg = 1, nimages, im, i, j, ii, jj, dx, dy, nx, ny;
	FILE	*ipf, *offsetf;
	int	N1, N2, *ioffset, *joffset, pixtype, drawbox, offsetsarg, useoffsetsfile, hasgrid;
	fitsheader	*fits;
	char	argstring[COM_LENGTH], errorstring[1024], line[1024];
	float	**fin, **fout, xoff, yoff;
	int	Nin1, Nin2, magicinit = 0;

	if (argc < 3)
		error_exit(usage);
	nimages = argc - 3;

	/* defaults */
	useoffsetsfile = 0;
	drawbox = 0;
	offsetsarg = 0;
	hasgrid = 0;

	/* process arguments */
        while (arg < argc) {
                if (argv[arg][0] == '-') {
                        switch (argv[arg++][1]) {
                                case 'b':
					drawbox = 1;
                                        nimages -= 1;
					break;
                                case 'M':
					magicinit = 1;
                                        nimages -= 1;
					break;
                                case 'f':
					offsetsarg = arg++;
                                        nimages -= 2;
					useoffsetsfile = 1;
					break;
                                case 'd':
                                        if (1 != sscanf(argv[arg++], "%d", &dx))
                                                error_exit(usage);
                                        if (1 != sscanf(argv[arg++], "%d", &dy))
                                                error_exit(usage);
					hasgrid = 1;
                                        nimages -= 3;
                                        break;
                                default:
                                        error_exit(usage);
                                        break;
                        }
                } else {
                        break;
                }
        }
	if ((argc - arg) < 3)
		error_exit(usage);

	if (useoffsetsfile) {
		sscanf(argv[arg++], "%d", &N1);
		sscanf(argv[arg++], "%d", &N2);
		ioffset = (int *) calloc(nimages, sizeof(int));
		joffset = (int *) calloc(nimages, sizeof(int));
		offsetf = fopen(argv[offsetsarg], "r");
		if (!offsetf) 
			error_exit("album: unable to open offsetsfile for input\n");
		fgets (line, 1024, offsetf);
		for (im = 0; im < nimages; im++) {
			if (!fgets (line, 1024, offsetf))
				error_exit("album: problem reading from offsets.out\n");
			if (2 != sscanf(line, "%f %f", &xoff, &yoff))
				error_exit("album: problem reading from offsets.out\n");
			joffset[im] = xoff;
			ioffset[im] = yoff;
		}
	} else {
		sscanf(argv[arg++], "%d", &nx);
		sscanf(argv[arg++], "%d", &ny);
		ioffset = (int *) calloc(nx * ny, sizeof(int));
		joffset = (int *) calloc(nx * ny, sizeof(int));
		if (!hasgrid) {
	        	ipf = fopen(argv[argc - nimages], "r");
                	if (!ipf) {
                        	sprintf(errorstring, "album: failed to open %s\n", 
                      	          argv[argc - nimages]);
                        	error_exit(errorstring);
               		}
			fits = readfitsheader(ipf);
			/* set_fits_ipf(ipf);
			read_fits_head(&dx, &dy, &comc, comv);*/
			fclose(ipf);
			dx = fits->n[0];
			dy = fits->n[1];
		}			
		N1 = nx * dx;
		N2 = ny * dy;
		for (i = 0; i < ny; i++) {
			for (j = 0; j < nx; j++) {
				ioffset[nx * i + j] = i * dy;
				joffset[nx * i + j] = j * dx;
			}
		}
	}


	allocFloatArray(&fout, N1, N2);
	if (magicinit) {
		for (i = 0; i < N2; i++) {
			for (j = 0; j < N1; j++) {
				fout[i][j] = FLOAT_MAGIC;
			}
		}
	}

	/* process files */
	for (im = 0; im < nimages; im++) {
	        ipf = fopen(argv[argc - nimages + im], "r");
                if (!ipf) {
                        sprintf(errorstring, "album: failed to open %s\n", 
                                argv[argc - nimages + im]);
                        error_exit(errorstring);
                }
		read2Dfloatimage(&fin, &Nin1, &Nin2, &fits, ipf);
		/*set_fits_ipf(ipf);
		fread_fits(&fin, &Nin1, &Nin2, &comc, comv);*/
		fclose(ipf);
		for (i = 0; i < Nin2; i++) {
			for (j = 0; j < Nin1; j++) {
				ii = ioffset[im] + i;
				jj = joffset[im] + j;
				if (ii >= 0 && ii < N2 && jj >= 0 && jj < N1) {
					if (fin[i][j] != FLOAT_MAGIC)
						fout[ii][jj] = fin[i][j];
					if (drawbox) {
						if (i == 0 || i == (Nin2 - 1) || j == 0 || j == (Nin1 - 1)) {
							fout[ii][jj] = SHRT_MAX;
						}
					}
				}
			}
		}
	}
	add_comment(argc, argv, fits);
	set2Dimagesize(fits, N1, N2);
	write2Dfloatimage(fout, fits);
	/* write_fits(fout, N1, N2, comc, comv); */
	exit(0);
}
