#define	usage "\n\n\n\
NAME\n\
	print_image -- prints a square image on the LaserWriter\n\
\n\
SYNOPSIS\n\
	print_image	fitsfile [option...] \n\
		-f fmin fmax	# max and min f values\n\
		-n		# don't print header info\n\
		-b		# draw a box oround image\n\
		-c comment	# use (quoted if multiword) comment as label\n\
		-p pw ph mx my	# set page width, height, margins (612, 792, 50, 50)\n\
		-P		# include setpagedevice\n\
		-g		# pipe output through gs to compress\n\
\n\
DESCRIPTION\n\
	\"print_image\" reads a fits file and sends\n\
	postscript file to stdout\n\
\n\
	The range of values may be specified with the -f option,\n\
	otherwise the range is 0 (=white) to 255 (=black).\n\
\n\
	Supports 3-color images.\n\
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


char		*psstring;
static	FILE	*opf;


void	print_im(int N1, int N2, int Ncolors, double (*gray)(), char *caption);	
void	setpagesize(int pagewidth, int pageheight, int pagemargin);
void	print_caption(char *caption);
void	psDrawChar(char aChar);
void	ps(char *string);
void	set_print_opf(FILE *thefile);

main(int argc, char *argv[])	
{
	int		i, j, arg = 2, N1, N2, Ncolors = 1, Naxes, N[3], docolor = 0, color;	
	int		dobox, icom;
	fitsheader	*fits;
	fitscomment	*com;
	char		*thecomment, nullcomment[2] = "";
	char		caption[1024], argstring[512];
	FILE		*imagef;
	int		ascii, noheader;
	int		caption_height, width, height;
	int		page_width, page_height, page_margin_x, page_margin_y, dosetpagedevice;
	int		psstringlen;
	float		*f, fmin, fmax, fval;
	unsigned char	***c;
		
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

	while (arg < argc) {
		if (*argv[arg] != '-')
			error_exit(usage);
		switch (*(argv[arg++]+1)) {
			case 'f':
				if (!sscanf(argv[arg++], "%f", &fmin))
					error_exit(usage);
				if (!sscanf(argv[arg++], "%f", &fmax))
					error_exit(usage);
				if (fmax == fmin)
					error_exit("print_image: fmin and fmax must be unequal\n");
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
				thecomment = argv[arg++];
				break;
			case 'p':
				if (!sscanf(argv[arg++], "%d", &page_width))
					error_exit(usage);
				if (!sscanf(argv[arg++], "%d", &page_height))
					error_exit(usage);
				if (!sscanf(argv[arg++], "%d", &page_margin_x))
					error_exit(usage);
				if (!sscanf(argv[arg++], "%d", &page_margin_y))
					error_exit(usage);
				break;
			case 'P':
				dosetpagedevice = 1;
				break;
			case 'g':
				opf = popen("gs -q -sDEVICE=pswrite -sOutputFile=- -dNOPAUSE -", "w");
				if (!opf) {
					error_exit("print_image: failed to open gs pipe for output\n");
				}
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	page_width -= 2 * page_margin_x;
	page_height -= 2 * page_margin_y;
	

	if (!strcmp(argv[1], "-u"))
		error_exit(usage);
	if (!strcmp(argv[1], "-"))
		imagef = stdin;
	else
		imagef = fopen(argv[1], "r");
	if (!imagef)
		error_exit("print_image: unable to open fits file for input\n");
	fits = readfitsheader(imagef);
	N1 = fits->n[0];
	N2 = fits->n[1];
	if (fits->ndim > 2) { 
		if (fits->n[2] != 3) {
			error_exit("print_image: I can only deal with 3-color images\n");
		} else {
			docolor == 1;
			Ncolors = 3;
		}
	}

	f = (float *) calloc(N1, sizeof(float));
	c = (unsigned char ***) calloc(Ncolors, sizeof(unsigned char **));
	for (color = 0; color < Ncolors; color++) {
		c[color] = (unsigned char **) calloc(N2, sizeof(unsigned char *));
		for (i = 0; i < N2; i++) {
			c[color][i] = (unsigned char *) calloc(N1, sizeof(unsigned char));
			readfitsline(f, fits);
			for (j = 0; j < N1; j++) {
				if (Ncolors != 1) {
					if (f[j] == FLOAT_MAGIC) {
						fval = 255;
					} else {
						fval = 256 * (f[j] - fmin) / (fmax - fmin);
					}
				} else {
					fval = 256 * (fmax - f[j]) / (fmax - fmin);
				}
				c[color][i][j] = (unsigned char) (fval >= 0 ? (fval < 256 ? fval : 255) : 0);
			}
		}
		if (dobox) {
			for (i = 0; i < N1; i++) {
				c[color][0][i] = c[color][N2 - 1][i] = (unsigned char) 0;
			}
			for (i = 0; i < N2; i++) {
				c[color][i][0] = c[color][i][N1 - 1] = (unsigned char) 0;
			}
		}
	}


	argsToString(argc, argv, argstring);
	sprintf(caption, "%s\n", argstring);
	if (noheader)
		strcpy(caption, thecomment);

	psstringlen = (2 * Ncolors * N1 + 2 > 1024 ? 2 * Ncolors * N1 + 2 : 1024);
	psstring = (char *) calloc(psstringlen, sizeof(char));

	if (strlen(caption)) {
		caption_height = 50;
	} else {
		caption_height = 0;
	}

	if (((float) N2 / (page_height - caption_height)) > ((float) N1 / page_width)) {
		height = page_height - caption_height;
		width = N1 * height / N2;
	} else {
		width = page_width;
		height = N2 * width / N1;
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

	fprintf(opf, "/picstr %d string def\n", N1 * Ncolors);
	ps("0 0 translate");
	
	fprintf(opf, "%d %d scale\n", (int) width, (int) height);
	switch (Ncolors) {
		case 1:
			sprintf(psstring, 
				"/drawImage {%d %d 8 [%d 0 0 %d 0 0] {currentfile picstr readhexstring pop} image} def", 
				N1, N2, N1, N2);
			break;
		case 3:
			sprintf(psstring, 
				"/drawImage {%d %d 8 [%d 0 0 %d 0 0] {currentfile picstr readhexstring pop} false 3 colorimage} def", 
				N1, N2, N1, N2);
			break;
		default:
			error_exit("print_im: I only know how to deal with one or three colors\n");
			break;
	}
	ps(psstring);
	ps("drawImage");
	for (i = 0; i < N2; i++) {
		for(j = 0; j < N1; j++) {
			for (color = 0; color < Ncolors; color++) {
				if (c[color][i][j] < 16) {
					fprintf(opf, "0%X", (int) c[color][i][j]);
				} else {
					fprintf(opf, "%X", (int) c[color][i][j]);
				}
			}
		}
		fprintf(opf, "\n");
	}
	ps("grestore");
	ps("showpage");

	exit(0);
}

/*
double	gray(int i, int j, int color)
{
	return((double) (f[0][i][j] - fmin) / (double) (fmax - fmin));
}

double	shade(int i, int j, int color)
{
	return((double) (fmax - f[color][i][j]) / (double) (fmax - fmin));
}

*/


void		print_im(int N1, int N2, int Ncolors, double (*gray)(), char *caption)	
{
	char		pixval[256][2];
	char		hexchar[16] = {'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};
	int		i, j, color, fval, caption_height;
	double		width, height;
	
	/* set up pixel values as hex digit pairs */
	for (i = 0; i < 16; i++) {
		for (j = 0; j < 16; j++) {
			pixval[255 - (16 * i + j)][0] = hexchar[i];
			pixval[255 - (16 * i + j)][1] = hexchar[j];
		}
	}	

	

	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {
			for (color = 0; color < Ncolors; color++) {
				fval = 256 * gray(i, j, color);	
				fval = (fval >= 0 ? (fval < 256 ? fval : 255) : 0);
				psstring[2 * (j * Ncolors + color)] = pixval[fval][0];
				psstring[2 * (j * Ncolors + color) + 1] = pixval[fval][1];
			}
		}
		psstring[2 * N1 * Ncolors] = '\0';
		ps(psstring);
	}
}

void	ps(char *string)
{
	fprintf(opf, "%s\n", string);
}


void	psDrawChar(char aChar)
{
	switch (aChar) {
		case '(':
			sprintf(psstring, "(\\50) show");
			break;
		case ')':
			sprintf(psstring, "(\\51) show");
			break;
		default:
			sprintf(psstring, "(%c) show", aChar);
			break;
	}
	ps(psstring);
}



void	print_caption(char *caption) {
	float	size = 10, space = 13, h = 13, v = -28;
	int		c, len;
	len = strlen(caption);
	
	ps("/Times-Roman findfont   10.000 scalefont setfont");
	sprintf(psstring, "%f %f moveto", h, v);
	ps(psstring);
	for (c = 0; c < len; c++)
		switch (caption[c]) {
			case '\n':
				v -= space;
				sprintf(psstring, "%f %f moveto", h, v);
				ps(psstring);
				break;
			default:
				psDrawChar(caption[c]);
		}
}



void	set_print_opf(FILE *thefile)
{
	opf = thefile;
}




