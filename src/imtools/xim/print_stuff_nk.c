/**
 ** print_stuff.c functions for printing images as postscript
 **
 **		public use:
 ** 		int		print_im(int N1, int N2, int Ncolors, double (*gray)(), char *caption);	
 **
 **		private use:
 **			void	print_caption(char *caption)
 **			void	psDrawChar(char aChar)
 **			void	ps(char *string)
 **
 **		print_im() is passed a function which returns gray(i,j,col) on range 0-1.0
 **		where col = 0,1,2 for r,g,b
 **/
 
 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <signal.h>

#include "print_stuff.h"
#include "../../utils/error.h"
#include "../../imlib/fits.h"

#define	PAGE_WIDTH			512
#define	PAGE_HEIGHT			692
#define	PAGE_MARGIN			50
#define	EPSF_A				50
#define	EPSF_B				50


static	int	page_width 	= PAGE_WIDTH;
static	int	page_height 	= PAGE_HEIGHT;
static	int	page_margin 	= PAGE_MARGIN;
static int	gdosetpagedevice = 0;


int			psstringlen;
char		*psstring;

static	FILE	*opf = stdout;

void		print_im(int N1, int N2, int Ncolors, double (*gray)(), char *caption)	
{
	char		pixval[256][2];
	char		hexchar[16] = {'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};
	int		i, j, color, fval, caption_height;
	char		boundingboxline[1024];
	double		width, height;
	
	/* set up pixel values as hex digit pairs */
	for (i = 0; i < 16; i++) {
		for (j = 0; j < 16; j++) {
			pixval[255 - (16 * i + j)][0] = hexchar[i];
			pixval[255 - (16 * i + j)][1] = hexchar[j];
		}
	}	
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


/*
	standard postscript header
	ps("\%!PS");
*/

/*	epsf header */
	ps("\%!PS-Adobe-2.0 EPSF-2.0");
	ps("\%\%Creator: imcat");
	sprintf(boundingboxline, "\%\%\%\%BoundingBox: %d %d %d %d", 
		page_margin, page_margin,
		(int) ceil(width + page_margin), 
		(int) ceil(height + page_margin + caption_height));
	ps(boundingboxline);
	ps("\%\%EndComments");

	/* add the PageSize directive */
	if (gdosetpagedevice) {
		sprintf(psstring, "<</PageSize [%d %d] >> setpagedevice", 
			page_width + 2 * page_margin, 
			page_height + 2 * page_margin);
		ps(psstring);
	}

	/* save the previous state of the printer */
	ps("gsave");

	/* set the origin at bot left of image square */
	sprintf(psstring, "%d %d translate", 
		page_margin, page_margin + caption_height);
	ps(psstring);

	print_caption(caption);

	sprintf(psstring, "/picstr %d string def", N1 * Ncolors);
	ps(psstring);
	ps("0 0 translate");
	
	sprintf(psstring, "%d %d scale", (int) width, (int) height);
	ps(psstring);
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
	ps("grestore");
	ps("showpage");
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



void	setpagesize(int pagewidth, int pageheight, int pagemargin, int dosetpagedevice)
{
	page_width = pagewidth;
	page_height = pageheight;
	page_margin = pagemargin;
	gdosetpagedevice = dosetpagedevice;
}

