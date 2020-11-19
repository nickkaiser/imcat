/*
#include 	<Xm/MainW.h>
#include 	<Xm/Label.h>
*/
#include 	<stdio.h>
#include	<stdlib.h>
#include	"imview_stuff.h"
#include	"../../imlib/fits.h"
#include	"../../utils/error.h"

#define		MAX_GRAYS	128
#define		SIZE_GOAL	512

/* global declarations */
static	Display 	*mydisplay;
static	int 		myscreen, depth;
static	GC		mygc;
static	Colormap	cmap;
static	Window		mywindow;
static	XColor		gray[MAX_GRAYS];
static	unsigned long	pixel[MAX_GRAYS];
static	int		n_grays;
extern	Widget		main_w;
extern	int		ascii;

void		alloc_grays(int *n_cells)
{
	int		i;
	XVisualInfo	info;
	Visual		*visual;
	int		class = PseudoColor;
	int		plane_masks[1];

	mydisplay = XtDisplay(main_w);
	myscreen = DefaultScreen(mydisplay);
	mywindow = DefaultRootWindow(mydisplay);	
	mygc = XCreateGC(mydisplay, mywindow, 0, 0);

	depth = DefaultDepth(mydisplay, myscreen);
	cmap = DefaultColormap(mydisplay, myscreen);
	visual = DefaultVisual(mydisplay, myscreen);
	if (depth == 1) {
		fprintf(stderr, "this program won't work on a monochrome display\n");
		exit(0);
	}

	if (!XMatchVisualInfo(mydisplay, myscreen, depth, class, &info)) {
		fprintf(stderr, "wrong display type\n");
		exit(0);
	}
	
	n_grays = MAX_GRAYS;
	while (1) {
		if (XAllocColorCells(mydisplay, cmap, False, plane_masks,
					0, pixel, n_grays))
			break;
		n_grays--;
	}	
	if (n_grays == 0) {
		fprintf(stderr, "couldn't allocate any grays\n");
		exit(0);
	}
/*	fprintf(stderr, "allocated %d color cells\n", n_grays);*/
	*n_cells = n_grays;
}




void	set_grays(int n_cells, int n_lo, int n_hi)
{
	int	i, inverse = 0, temp;
	long	f, df;

	if (n_hi == n_lo)
		if (n_hi > 0)
			n_lo = n_hi - 1;
		else
			n_hi = n_lo + 1;
	if (n_hi < n_lo) {
		inverse = 1;
		temp = n_hi; n_hi = n_lo; n_lo = temp;
	}

	df = 65535 / (n_hi - n_lo);

	for (i = 0; i < n_cells; i++)
		{
			if (i < n_lo)
				f = (inverse ? 65535 : 0);
			else
				if (i < n_hi)
					f = (inverse ?
						65535 - (i - n_lo) * df :
						(i - n_lo) * df);
				else
					f = (inverse ? 0 : 65535);
			gray[i].pixel = pixel[i];
			gray[i].red = gray[i].green = gray[i].blue = f;
			gray[i].flags = DoRed | DoGreen | DoBlue;				
		}
	XStoreColors(mydisplay, cmap, gray, n_cells);
}



Pixmap	make_pixmap(int fMin, int fMax, int n_cells) 
{
	unsigned int	imWidth, imHeight, format; 
	int		i, j, ii, jj, bitmap_pad, color_index;	
	Pixmap		thepixmap;
	XImage 		*image;
	int		N, comc, scale = 1;
	char		*comv[MAX_COMMENTS];
	short		*f;

	read_fits_head(&N, &N, &comc, comv);
	if (N < SIZE_GOAL)
		scale = SIZE_GOAL / N;

	imWidth = imHeight = N * scale;
	f = (short *) calloc(N, sizeof(short));
	if (!f)
		error_exit("make_pixmap: memory allocation error\n");

	if (depth == 1)
		{
			format = XYPixmap;
			bitmap_pad = 32;
		}
	else
		{
			format = ZPixmap;
			bitmap_pad = 8;
			if (depth > 8) bitmap_pad = 32;
		}
	image = XCreateImage(mydisplay,
			DefaultVisual(mydisplay, myscreen),
			depth, format, 0, 0, imWidth, imHeight,
			bitmap_pad, 0);
	if (image == 0)
		{
			printf("image allocation failure\n"); 
			exit(1);
		}
	image->data = (char *) malloc(image->bytes_per_line * imHeight);
	if (image->data == 0)
		{
			printf("image memory allocation failure\n"); 
			exit(1);
		}

	for (j = N - 1; j >= 0; j--) {
		read_fits_line(f, N);
		for (i = 0; i < N; i++)
			{
				color_index = n_cells * (f[i] - fMin)
					/ (float) (fMax - fMin);
				if (color_index >= n_cells) color_index = n_cells - 1;
				if (color_index <  0) color_index = 0;
				if (scale == 1)
					XPutPixel(image, i, j, gray[color_index]);
				else
					for (ii = scale * i; ii < scale * (i + 1); ii++)
						for (jj = scale * j; jj < scale * (j + 1); jj++)
							XPutPixel(image, ii, jj, gray[color_index]);
			}
	}
/*	printf("image painting done....\n");*/
	free(f);

	/* pixmap creation */
	thepixmap = XCreatePixmap(mydisplay, mywindow, imWidth, imHeight, depth);
/*	printf("pixmap created....\n");*/

	XPutImage(mydisplay,thepixmap, mygc, 
			image, 0, 0, 0, 0, imWidth, imHeight);
	return (thepixmap);
}



