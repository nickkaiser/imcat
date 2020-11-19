/*
 * fitstops.c - main program
 */


#define	usage "\n\n\n\
NAME\n\
	fitstops -- convert image from fits to postscript\n\
\n\
SYNOPSIS\n\
	fitstops fitsfile [option...] \n\
		-f fmin fmax	# max and min f values\n\
		-n		# don't print header info\n\
		-b		# draw a box oround image\n\
		-c comment	# use (quoted if multiword) comment as label\n\
		-p pw ph mx my	# set page width, height, xmargin, ymargin (612, 792, 50, 50)\n\
		-P		# include setpagedevice command\n\
		-g		# pipe output through gs to compress\n\
		-C nc lw	# draw nc contours of width lw\n\
		-s		# draw 3D hidden line surface plot (alt = az = 30.0; zfac = 3.0; zoff = -0.2)\n\
		-S alt az zfac zoff	# draw 3D hidden line surface plot\n\
\n\
DESCRIPTION\n\
	\"fitstops\" reads fitsfile and, by default, sends\n\
	postscript gray scale or color image to stdout\n\
\n\
	The range of values may be specified with the -f option,\n\
	otherwise the range is 0 (=white) to 255 (=black).\n\
\n\
	Supports 1 and 3-color images.\n\
\n\
	If fitsfile = '-' then we read from stdin.\n\
\n\
	Default page size info is 612x792 = (8.5x11)in with 50 pt margin\n\
	so the actual inked area is 512x692. Use -p option to change this.\n\
\n\
	The -P option is provided to include a 'setpagedevice' command giving\n\
	the physical total page size.  This is used for big prints on\n\
	the designjet, but seems to be problematic with latex epsf handling.\n\
	Do not use this with -g option.\n\
\n\
	With the -g option we pipe the output through\n\
		gs -q -sDEVICE=pswrite -sOutputFile=- -dNOPAUSE -\n\
	which will result in a much smaller output file.\n\
	Do not use this with -P option.\n\
\n\
	Use the -s, -S options to generate a 3D hidden line surface plot using\n\
	Tonry's mongo routine.  The simple -s option uses default parameters\n\
	which work reasonably well for a positive peak of height unity near the\n\
	origin.  Use the -S option to fiddle with the parameters alt, az, zfac, zoff\n\
	where:\n\
		alt, az  = altitude azimuth viewing angle in degrees.\n\
		zfac = scaling of z-axis (roughly 3.0 / data max)\n\
		zoff = offset of z-origin.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"		

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../utils/args.h"
#include "psutils.h"
#include "makec.h"
#include "printimage.h"
#include "contourplot.h"
#include "surfaceplot.h"
#include "tonry3d.h"

main(int argc, char *argv[])	
{
	int		i, j, arg = 2, N1, N2, Ncolors = 1, Naxes, N[3], docolor = 0, color;	
	int		dobox, icom;
	fitsheader	*fits;
	fitscomment	*com;
	char		*thecomment, nullcomment[2] = "";
	char		caption[1024], argstring[512], *srcfilename, *flag;
	FILE		*imagef;
	int		ascii, noheader, docontours, icont, ncont;
	int		caption_height, width, height;
	int		page_width, page_height, page_margin_x, page_margin_y, dosetpagedevice;
	int		psstringlen;
	int		dosurface;
	float		*f, fmin, fmax, fval, *fcont, lw;
	unsigned char	***c;
	FILE		*opf;
	float		alt, az, zfac, zoff;
		
	/* defaults */
	fmin = 0.0;
	fmax = 255.0;
	dobox = 0;
	thecomment = nullcomment;
	page_width = 612;
	page_height = 792;
	page_margin_x = 50;
	page_margin_y = 50;
	dosetpagedevice = 0;
	ascii = 0;
	noheader = 0;
	opf = stdout;
	docontours = 0;
	dosurface = 0;
	alt = 30.0;
	az = 30.0;
	zfac = 3.0;
	zoff = -0.2;

	argsinit(argc, argv, usage);
	srcfilename = getargs();
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'f':
				fmin = getargf();
				fmax = getargf();
				if (fmax == fmin)
					error_exit("fitstops: fmin and fmax must be unequal\n");
				break;
			case 'a':
				ascii = 1;
				break;
			case 'n':
				noheader = 1;
				break;
			case 'b':
				dobox = 1;
				break;
			case 'c':
				noheader = 1;
				thecomment = getargs();
				break;
			case 'p':
				page_width  = getargi();
				page_height = getargi();
				page_margin_x = getargi();
				page_margin_y = getargi();
				break;
			case 'P':
				dosetpagedevice = 1;
				break;
			case 'g':
				opf = popen("gs -q -sDEVICE=pswrite -sOutputFile=- -dNOPAUSE -", "w");
				if (!opf) {
					error_exit("fitstops: failed to open gs pipe for output\n");
				}
				break;
			case 'C':
				docontours = 1;
				ncont = getargi();
				lw = getargf();
				break;
			case 'S':
				alt = getargf();
				az  = getargf();
				zfac = getargf();
				zoff = getargf();
			case 's':
				dosurface = 1;
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	set_print_opf(opf);

	page_width -= 2 * page_margin_x;
	page_height -= 2 * page_margin_y;
	

	if (!strcmp(argv[1], "-u"))
		error_exit(usage);
	if (!strcmp(argv[1], "-"))
		imagef = stdin;
	else
		imagef = fopen(argv[1], "r");
	if (!imagef)
		error_exit("fitstops: unable to open fits file for input\n");
	fits = readfitsheader(imagef);
	N1 = fits->n[0];
	N2 = fits->n[1];
	if (fits->ndim > 2) { 
		if (fits->n[2] != 3) {
			error_exit("fitstops: I can only deal with 3-color images\n");
		} else {
			docolor == 1;
			Ncolors = 3;
		}
	}



	argsToString(argc, argv, argstring);
	sprintf(caption, "%s\n", argstring);
	if (noheader)
		strcpy(caption, thecomment);


	if (strlen(caption)) {
		caption_height = 50;
	} else {
		caption_height = 0;
	}

	if (dosurface) {
		height = page_height - caption_height;
		width = page_width;
	} else {
		if (((float) N2 / (page_height - caption_height)) > ((float) N1 / page_width)) {
			height = page_height - caption_height;
			width = N1 * height / N2;
		} else {
			width = page_width;
			height = N2 * width / N1;
		}
	}


	ps("\%!PS-Adobe-2.0 EPSF-2.0");
	ps("\%\%Creator: imcat");

	fprintf(opf, "\%\%\%\%BoundingBox: %d %d %d %d\n", 
		page_margin_x, page_margin_y,
		(int) ceil(width + page_margin_x), 
		(int) ceil(height + page_margin_y + caption_height));
	ps("\%\%EndComments");

	/* add the PageSize directive */
	if (dosetpagedevice) {
		fprintf(opf, "<</PageSize [%d %d] >> setpagedevice\n", 
			page_width + 2 * page_margin_x, 
			page_height + 2 * page_margin_x);
	}
	/* save the previous state of the printer */
	ps("gsave");

	/* set the origin at bot left of image square */
	fprintf(opf, "%d %d translate\n", 
		page_margin_x, page_margin_y + caption_height);

	print_caption(caption);

	if (dosurface) {
		fprintf(opf, "%10.3f setlinewidth\n", 32.0 / N1);
		fprintf(opf, "/n {newpath} bind def\n");
		fprintf(opf, "/m {moveto} bind def\n");
		fprintf(opf, "/l {lineto stroke} bind def\n");
		surfaceplot(fits, (float) width, (float) height, alt, az, zfac, zoff);
	} else if (docontours) {
		fprintf(opf, "%10.3f setlinewidth\n", lw);
		fprintf(opf, "/n {newpath} bind def\n");
		fprintf(opf, "/m {moveto} bind def\n");
		fprintf(opf, "/l {lineto stroke} bind def\n");
		fcont = (float *) calloc(ncont, sizeof(float));
		for (icont = 0; icont < ncont; icont++) {
			fcont[icont] = fmin + icont * (fmax - fmin) / ncont;
		}
		contourplot(fits, N1, N2, ncont, fcont, (float) width, (float) height);
	} else {
		c = makecarray(fits, N1, N2, Ncolors, fmin, fmax, dobox);
		printimage(Ncolors, N1, N2, c, width, height);
	}

	ps("grestore");
	ps("showpage");

	if (opf == stdout) {
		exit(0);
	} else {
		exit(pclose(opf));
	}
}








