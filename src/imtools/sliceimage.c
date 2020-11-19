#define usage "\n\n\n\
NAME\n\
	sliceimage - cut a 2-D FITS image into a grid of sub-sections\n\
\n\
SYNOPSIS\n\
	sliceimage [-u] [-i] [-l] [-f|-c format] nx ny\n\
\n\
DESCRIPTION\n\
	In its default mode, sliceimage reads a 2D image fin[N2][N1], whose fast\n\
	dimension N1 must be a multiple of nx and whose slow index\n\
	N2 must be a multiple of ny, from stdin and sends to stdout\n\
	a 4-D image fout[ny][nx][N2/ny][N1/nx] such that\n\
\n\
		fout[y][x][Y][X] = fin[y * N2 / ny + Y][x + N1 / nx + X]\t\n\
\n\
	This slices a single image into a grid of contiguous patches.\n\
\n\
	With -f option, we write the results to a set of files with names\n\
	generated according to the format string 'format'.\n\
\n\
	With -i option it performs the inverse operation. The arguments\n\
	nx and ny are ignored.\n\
\n\
	With -l option it outputs instead a 3-dimensional image\n\
\n\
		fout[y * nx + x][Y][X] = fin[y * N2 / ny + Y][x + N1 / nx + X]\t\n\
\n\
	and writes nx and ny to the FITS header as SLICE_NX and SLICE_NY.\n\
\n\
	With -f option, we write the results to a set of files with names\n\
	generated according to the format string 'format'.  With -l option\n\
	the output filenames are generated as\n\
\n\
		l = y * nx + x;\n\
		sprintf(filename, format, l, l, l, l, l, l);\n\
\n\
	so you can use the index l multiple times in the filename (or command)\n\
	otherwise the filenames are generated by\n\
\n\
		sprintf(filename, format, y, x);\n\
\n\
	With -c option, the sub-images are piped through a command which is\n\
	generated from the string format just as with the -f option.\n\
\n\
	With -u option it sends this man page to stderr and exits\n\
	with abnormal status.\n\
\n\
BUGS\n\
	You cannot use the -c option to send the output of the subprocesses\n\
	sequentially to stdout.  You have to redirect the output of these\n\
	processes to files.  However, you can allways collect the output from\n\
	these files with stackplanes, as, for example:\n\
\n\
		sliceimage -l -c 'foo > tmp.%%d.fits' 3 3 ; stackplanes tmp.[0-8].fits\n\
\n\
SEE ALSO\n\
	flatten(1), getplane(1), getplanes(1), stackplanes(1).\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/args.h"

int		main(int argc, char *argv[])	
{
	fitsheader	*fits;
	fitscomment *thecomment;
	char		*f, **finv;
	int		inverseopmode, fits3Dmode, writefilesmode, applycommandmode;
	char		*flag, filename[1024], *formatstring;
	int		NX, NY, nx, ny, dx, dy, ix, iy, l, bytesperpixel;
	FILE		*opf = stdout;

	/* defaults */
	writefilesmode = 0;
	inverseopmode = 0;
	fits3Dmode = 0;
	applycommandmode = 0;

	/* parse args */
	argsinit(argc, argv, usage);
	while (nextargtype() == FLAG_ARG) {
		flag = getflag();
		switch (flag[0]) {
			case 'i':
				inverseopmode = 1;
				break;
			case 'l':
				fits3Dmode = 1;
				break;
			case 'f':
				writefilesmode = 1;
				formatstring = getargs();
				break;
			case 'c':
				applycommandmode = 1;
				formatstring = getargs();
				break;
			default:
				error_exit(usage);
		}
	}
	if (!inverseopmode) {
		nx = getargi();
		ny = getargi();
	}
	if (inverseopmode && (writefilesmode || applycommandmode)) {
		error_exit("sliceimage : you can't use -f (or -c) and -i options together");
	}
	if (writefilesmode && applycommandmode) {
		error_exit("sliceimage : you can't use -f and -c options together");
	}

	/* read header and add hisotry */
	fits = readfitsheader(stdin);
	add_comment(argc, argv, fits);

	/* write the array dimensions if we're making a 3D FITS file */
	if (fits3Dmode && !inverseopmode) {
		removenamedcomments("SLICE_NX", fits);
		appendcomment(newnumericcomment("SLICE_NX", (double) nx, "sliceimage nx argument"), fits);
		removenamedcomments("SLICE_NY", fits);
		appendcomment(newnumericcomment("SLICE_NY", (double) ny, "sliceimage ny argument"), fits);
	}

	/* get the pixel size from the header */
	bytesperpixel = pixsize(fits->extpixtype);

	/* we're just passing the data through, so */
	fits->opbyteorder = fits->ipbyteorder;

	/* do the biz */
	if (!inverseopmode) {
		if (fits->ndim != 2) {
			error_exit("sliceimage: input image must be 2-dimensional\n");
		}
		if (fits->n[0] % nx || fits->n[1] %ny) {
			error_exit("sliceimage: input image size must be exact multiple of ny, nx\n");
		}
		NX = fits->n[0];
		NY = fits->n[1];
		dx = fits->n[0] /= nx;
		dy = fits->n[1] /= ny;
		if (!(writefilesmode || applycommandmode)) {
			if (fits3Dmode) {
				fits->ndim = 3;
				fits->n[2] = nx * ny;
			} else {
				fits->ndim = 4;
				fits->n[2] = nx;
				fits->n[3] = ny;
			}
			writefitsheader(fits);
		}
		f = (char *) calloc(NX * dy * bytesperpixel, sizeof(char));
		for (iy = 0; iy < ny; iy++) {
			fread(f, sizeof(char), NX * dy * bytesperpixel, stdin);
			for (ix = 0; ix < nx; ix++) {
				if (writefilesmode || applycommandmode) {
					if (fits3Dmode) {
						l = iy * nx + ix;
						sprintf(filename, formatstring, l, l, l, l, l, l);
					} else {
						sprintf(filename, formatstring, iy, ix);
					}
					fits->opstream = opf = (writefilesmode ? fopen(filename, "w") : popen(filename, "w"));
					writefitsheader(fits);
				}
				for (l = 0; l < dy; l++) {
					fwrite(f + (l * NX + ix * dx) * bytesperpixel, dx * bytesperpixel, sizeof(char), opf);
				}
				if (writefilesmode || applycommandmode) {
					writefitstail(fits);
					(writefilesmode ? fclose(opf) : pclose(opf));
				}
			}
		}
		free(f);
	} else {
		switch (fits->ndim) {
			case 4:
				ny = fits->n[3];
				nx = fits->n[2];
				break;
			case 3:
				nx = (int) getnumericvalue(getcommentbyname("SLICE_NX", fits));
				ny = (int) getnumericvalue(getcommentbyname("SLICE_NY", fits));
				if (fits->n[2] != (nx * ny)) {
					error_exit("sliceimage: bad SLICE_N[XY] header values\n");
				}
				break;
			default:
				error_exit("sliceimage: input image must be 3 or 4-dimensional\n");
		}
		dx = fits->n[0];
		dy = fits->n[1];
		fits->n[0] *= nx;
		fits->n[1] *= ny;
		fits->ndim = 2;		
		writefitsheader(fits);
		finv = (char **) calloc(nx, sizeof(char *));
		for (ix = 0; ix < nx; ix++) {
			finv[ix] = (char *) calloc(dx * dy * bytesperpixel, sizeof(char));
		}
		for (iy = 0; iy < ny; iy++) {
			for (ix = 0; ix < nx; ix++) {
				fread(finv[ix], sizeof(char), dx * dy * bytesperpixel, stdin);
			}
			for (l = 0; l < dy; l++) {
				for (ix = 0; ix < nx; ix++) {
					fwrite(finv[ix] + (l * dx) * bytesperpixel, sizeof(char), dx * bytesperpixel, stdout);
				}
			}
		}
		for (ix = 0; ix < nx; ix++) {
			free(finv[ix]);
		}
		free(finv);	
	}
	if (!writefilesmode && !applycommandmode) {
		writefitstail(fits);
	}


	exit(0);
}
