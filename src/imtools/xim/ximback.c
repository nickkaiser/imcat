/*
 * ximback.c
 *
 */

#include 	<X11/Xlib.h>
#include 	<X11/Xutil.h>
#include 	<X11/Intrinsic.h>
#include 	<X11/StringDefs.h>
#include 	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<limits.h>
#include	"ximfront.h"
#include	"ximback.h"
#include	"ximcolor.h"
#include	"../../imlib/fits.h"
#include	"../../utils/error.h"
#include	"../../utils/stats_stuff.h"

#define		MAX_SHADES	4096

#define	MAX(x,y) (((x) > (y)) ? (x) : (y))
#define	MIN(x,y) (((x) < (y)) ? (x) : (y))

#define MAGIC SHORT_MAGIC

/* local globals */
static	int	level, nlevels, ncolors;
static	int	lastuseri, lastuserj, blocksize = 1;
static	int 		myscreen, depth;
static	GC		mygc;
static	Window		mywindow;
static	XImage		*theimage;
static	int		imwidth, imheight;
static	int		zoomi, zoomj;
static	Pixmap		thezoompixmap, pixmap[MAX_LEVELS];
char		caption[2048];

/* needed by ximcolor.c */
Colormap	cmap;
Display 	*mydisplay;
XColor		shade[MAX_SHADES];
unsigned 	long	pixel[MAX_SHADES];
unsigned int		n_shades;
short		**f[MAX_FITS_DIM][MAX_LEVELS];
int		 N1[MAX_LEVELS], N2[MAX_LEVELS];
fitsheader	*gfits;



int	xbgetvalue(int thename)
{
	switch (thename) {
		case LEVEL:
			return(level);
			break;
		case NLEVELS:
			return(nlevels);
			break;
		case NCOLORS:
			return(ncolors);
			break;
		case LASTUSERI:
			return(lastuseri);
			break;
		case LASTUSERJ:
			return(lastuserj);
			break;
		default:
			error_exit("xbgetvalue: bad name\n");
			break;
	}
}




void	xbsetvalue(int thename, int thevalue)
{
	switch (thename) {
		case LEVEL:
			level = thevalue;
			break;
		case LASTUSERI:
			lastuseri = thevalue;
			break;
		case LASTUSERJ:
			lastuserj = thevalue;
			break;
		default:
			error_exit("xbsetvalue: bad name\n");
			break;
	}
}



void		alloc_shades(void)
{
	int		i, color, shadesgoal;
	XVisualInfo	info;
	Visual		*visual;
	int		class = PseudoColor;
	unsigned long		plane_masks[1];

	mydisplay = XtDisplay(xfgetwidget(TOPLEVEL));
	myscreen = DefaultScreen(mydisplay);
	mywindow = DefaultRootWindow(mydisplay);	
	mygc = XCreateGC(mydisplay, mywindow, 0, 0);

	depth = DefaultDepth(mydisplay, myscreen);
	fprintf(stderr, "%d bit display\n", depth);
	shadesgoal = floor(0.5 + pow(2.0, (double) depth));
	fprintf(stderr, "shadesgoal = %d\n", shadesgoal);
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

	n_shades = shadesgoal;
	while (!XAllocColorCells(mydisplay, cmap, False, plane_masks, 0, pixel, n_shades))
		n_shades--;
	fprintf(stderr, "allocated %d color cells\n", n_shades);
	
	/* set up the color space */
	setcolorscheme();
}




void	readdataheader()
{
	int		level = 0, Naxes, *Naxis, nwords, color;
	char		word[7][64];
	int		goodsize;
	fitsheader	*fits;
	fitscomment	*com;
	
	fprintf(stderr,  "reading data header...\n");
	fits = readfitsheader(stdin);
	fits->intpixtype = SHORT_PIXTYPE;
	gfits = fits;
	Naxes = fits->ndim;
	Naxis = fits->n;

	com = fits->basecomment;
	while(com) {
                strncat(caption, com->name, NAME_LENGTH);
                strncat(caption, "= ", 2);
                strncat(caption, com->value, VALUE_LENGTH);
                strncat(caption, "\n", 1);
		com = com->next;
break;
        }

	N1[0] = Naxis[0];
	N2[0] = Naxis[1];
	ncolors = 1;
	if (Naxes == 3)
		ncolors = Naxis[2];
	if (Naxes > 3)
		error_exit("readdataheader: too many dimensions\n");
	fprintf(stderr, "ncolors = %d\n", ncolors);
	fprintf(stderr, "N1 = %d : N2 = %d\ncomments:\n", N1[0], N2[0]);
	com = fits->basecomment;
	while(com) {
		fprintf(stderr, "%s= %s\n", com->name, com->value);
		com = com->next;
        }
	goodsize = xfgetvalue(GOODSIZE);
	level = 0;
	while (N1[level] > goodsize || N2[level] > goodsize) {
		level++;
		N1[level] = N1[level - 1]  / 2;		/* round down on division */
		N2[level] = N2[level - 1]  / 2;
	}
	nlevels = level + 1;
	/* if image size < goodsize / 2 we want a unscrunched base level view */
	if (level == 0 && 2 * N1[0] <= goodsize && 2 * N2[0] <= goodsize) {
		blocksize = 1;
		while (2 * blocksize * N1[0] <= goodsize && 2 * blocksize * N2[0] <= goodsize)
			blocksize *= 2;
	}
	if (nlevels > MAX_LEVELS)
		error_exit("readdataheader: to many levels\n");
}


void	allocatedata(void)
{
	int	level, color;

	fprintf(stderr,  "allocating space for data...\n");
	for (color = 0; color < ncolors; color++)
		for (level = 0; level < nlevels; level ++)
			f[color][level] = alloc_f(N1[level], N2[level]);	
}


void	readdata(void)
{
	int		i, color;
	
	for (color = 0; color < ncolors; color++)
		for (i = 0; i < N2[0]; i++)
			readfitsline(f[color][0][i], gfits);
}


void	makescrunchedviews(void)
{
	int	i, j, ii, jj, level, color;
	double	fsum, nsum;

	fprintf(stderr,  "scrunching data...\n");
	for (color = 0; color < ncolors; color++) {
		for (level = 1; level < nlevels; level++) {
			N1[level] = N1[level - 1]  / 2;		/* round down on division */
			N2[level] = N2[level - 1]  / 2;
			for (i = 0; i < N2[level]; i++) {
				for (j = 0; j < N1[level]; j++) {
					fsum = nsum = 0;
					for (ii = 0; ii < 2; ii++) {
						for (jj = 0; jj < 2; jj++) {
							if (f[color][level - 1][2 * i + ii][2 * j + jj] != MAGIC) {
								nsum += 1.0;
								fsum += (float) f[color][level - 1][2 * i + ii][2 * j + jj];
							}
						}
					}
					if (nsum > 0)
						f[color][level][i][j] = (short) floor(0.5 + fsum / nsum);
					else
						f[color][level][i][j] = MAGIC;
				}
			}
		}
	}
}


/*
 * creat_pixmap():  creates a pixmap
 */
Pixmap	create_pixmap(int level) 
{
	Pixmap	thepixmap;

	thepixmap = XCreatePixmap(mydisplay, mywindow, blocksize * N1[level], blocksize * N2[level], depth);
	fprintf(stderr, "memory reserved for %4d x %4d x %2d pixmap....\n", 
		blocksize * N1[level], blocksize * N2[level], depth);
	return (thepixmap);
}

/*
 * fill_pixmap():  fills a pixmap
 */
void	fill_pixmap(void) 
{
	int		x, y, m, n, ii, jj, xx, yy, width, height;	
	unsigned long theshade;

	fprintf(stderr, "filling level %d pixmap....\n", level);
	
	width = imwidth / blocksize;
	height = imheight / blocksize;

	for (m = 0; m < N2[level] / height; m++) {
		for (n = 0; n < N1[level] / width; n++) {
			for (y = 0; y < height; y++) {
				ii = y + m * height;
				for (x = 0; x < width; x++) {
					jj = x + n * width;
					/* here is where we flip the y-axis: if you change this */
					/* be sure also to change GetPixel() */
					theshade = shade[color_index(level, N2[level] - ii - 1, jj)].pixel;
					for (xx = x * blocksize; xx < (x + 1) * blocksize; xx++) {
						for (yy = y * blocksize; yy < (y + 1) * blocksize; yy++) {
							XPutPixel(theimage, xx, yy, theshade);
						}
					}
				}
			}
			XPutImage(mydisplay, pixmap[level], mygc, 
				theimage, 0, 0, n * imwidth, m * imheight, imwidth, imheight);
		}
	}
}



void	makepreview(void) 
{
	int		i, j, step, color;

	fprintf(stderr, "making preview pixmap....\n");
	
	/* calculate the sampling rate */
	step = 1;
	for (i = 0 ; i < level; i++)
		step *= 2;
		
	for (color = 0; color < ncolors; color++) {
		for (i = 0; i < N2[level]; i++) {
			for (j = 0; j < N1[level]; j++) {
				f[color][level][i][j] = f[color][0][step * i][step * j];
			}
		}
	}
	fill_pixmap();
}



/*
 * makes a small XImage object which we will use to transfer
 * data to the pixmaps
 */
void	makeXImage(int N1, int N2)
{
	unsigned int	format; 
	int		bitmap_pad;	
	int	width, height;
	
	imwidth = width = N1 * blocksize;
	imheight = height = N2 * blocksize;

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
	theimage = XCreateImage(mydisplay,
			DefaultVisual(mydisplay, myscreen),
			depth, format, 0, 0, width, height,
			bitmap_pad, 0);
	if (theimage == 0)
		{
			printf("image allocation failure\n"); 
			exit(1);
		}
	theimage->data = (char *) malloc(theimage->bytes_per_line * height);
	if (theimage->data == 0)
		{
			printf("image memory allocation failure\n"); 
			exit(1);
		}
}


short	**alloc_f(int N1, int N2)
{
	int	i;
	short **f;
	
	f = (short **) calloc(N2, sizeof(short *));
	for (i = 0; i < N2; i++)
		f[i] = (short *) calloc(N1, sizeof(short));
	return (f);
}



void	GetPixelValue(int pixx, int pixy, int height, int width,
	int *useri, int *userj, short *pixval)
{
	int 	i, j, lev = level, color, step = 1, xsize, ysize;
	
	/* first we have to subtract half the border from the pixvalues */
	if (level >= 0) {
		xsize = blocksize * N1[level];
		ysize = blocksize * N2[level];
	} else {
		xsize = imwidth;
		ysize = imheight;
	}
	pixx -= (width - xsize) / 2;
	pixy -= (height - ysize) / 2;
	
	/* now divide by block size to get into user pixel coords */
	pixx /= blocksize;
	pixy /= blocksize;
	
	if (level >= 0) {
		i = pixy;
		j = pixx;
		while (lev > 0) {
			i *= 2;
			j *= 2;
			lev--;
		}
		if (pixx >= 0 && pixx < N1[level] && pixy >= 0 && pixy < N2[level]) {
			for (color = 0; color < ncolors; color++)
				pixval[color] = f[color][level][N2[level] - pixy - 1][pixx];
		} else {
			for (color = 0; color < ncolors; color++)
				pixval[color] = 0;
		}
		*useri = N2[0] - i - 1;
		*userj = j;
	} else {
		while (lev < 0) {
			step *= 2;
			lev++;
		}
		i = zoomi - (int) floor((double) (pixy - imheight / (2 * blocksize)) / (double) step);
		j = zoomj + (int) floor((double) (pixx - imwidth / (2 * blocksize)) /  (double) step);
		if (i >= 0 && i < N2[0] && j >= 0 && j < N1[0]) {
			for (color = 0; color < ncolors; color++)
				pixval[color] = f[color][0][i][j];	
		} else {
			for (color = 0; color < ncolors; color++)
				pixval[color] = 0;	
		}
		*useri = i;
		*userj = j;
	}
}




Pixmap	createzoompixmap(void)
{
	thezoompixmap = XCreatePixmap(mydisplay, mywindow, imwidth, imheight, depth);
}



void	fillzoompixmap(void)
{
	int	lev, step, x, y, i, j, xx, yy, width, height, izoom, jzoom;
	unsigned long	theshade;
	
	width = imwidth / blocksize;
	height = imheight / blocksize;
	
	/* set central point used by GetPixelValue() */
	zoomi = izoom = xbgetvalue(LASTUSERI);
	zoomj = jzoom = xbgetvalue(LASTUSERJ);
	
	step = 1;
	lev = level;
	while (lev < 0) {
		step *= 2;
		lev++;
	}
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			j = jzoom + (int) floor((double) (x - width / 2) / (double) step);
			i = izoom - (int) floor((double) (y - height / 2) / (double) step);
			if (i >= 0 && i < N2[0] && j >= 0 && j < N1[0]) {
				theshade = shade[color_index(0, i, j)].pixel;
			} else {
				theshade = shade[0].pixel;
			}
			for (xx = x * blocksize; xx < (x + 1) * blocksize; xx++) {
				for (yy = y * blocksize; yy < (y + 1) * blocksize; yy++) {
					XPutPixel(theimage, xx, yy, theshade);
				}
			}
		}
	}
	XPutImage(mydisplay, thezoompixmap, mygc, theimage, 0, 0, 0, 0, imwidth, imheight);
}

void	destroyzoompixmap(void)
{
	XFreePixmap(mydisplay, thezoompixmap);
}


int		myworkproc1(XtPointer client_data)
{
	int	color, slide;
	char label[20];

	/* readdataheader() sets up globals nlevels, blocksize, N1[], N2[], gfits*/
	readdataheader();
	if (ncolors > MAX_COLORS)
		error_exit("myworkproc1: too many colors\n");					 

	switch (ncolors) {
		case 1:
			setsliderval(0, 0.0);
			setsliderval(1, 1.0);
			XtUnmanageChild(xfgetwidget(SLIDER2));
			break;
		case 2:
			setsliderval(0, 1.0);
			setsliderval(1, 1.0 / A_MAX);
			setsliderval(2, 0.5);
			break;
		case 3:
			for (slide = 0; slide < 3; slide++) {
				setsliderval(slide, 1.0);
			}
			break;
		default:
			error_exit("myworkproc1: bad ncolours\n");
			break;
	}


	allocatedata();						/* now GetPixelValue will work ok */
	alloc_shades();			/* sets many globals like mydisplay etc. */
	for (slide = 0; slide < 3; slide++)			/* in case we print before wiggling sliders */
		newslidervalue(slide);
	set_shades();
	level = nlevels - 1;				/* most scrunched level */
	makeXImage(N1[level], N2[level]);	/* make an XImage the size of smallest view */

	/* now create the pixmaps */
	while (level >= 0) {
		pixmap[level] = create_pixmap(level);
		level --;
	}
	level = nlevels - 1;
	readdata();
	if (nlevels > 1) {
		makepreview();
		showpixmap();
	}
	return(1);
}



void    showpixmap(void)
{
        Cardinal                n;
        Arg             arg[5];
        
        n = 0;
	if (level >= 0) {
        	XtSetArg(arg[n], XtNheight, N2[level]);  n++;
        	XtSetArg(arg[n], XtNwidth, N1[level]); n++;
        	XtSetArg(arg[n], XtNbitmap, pixmap[level]); n++;
	} else {
        	XtSetArg(arg[n], XtNheight, imheight);  n++;
        	XtSetArg(arg[n], XtNwidth, imwidth); n++;
        	XtSetArg(arg[n], XtNbitmap, thezoompixmap); n++;
	}
        XtSetValues(xfgetwidget(PICWIDGET), arg, n);
}





int		myworkproc2(XtPointer client_data)
{
	makescrunchedviews();
	return(1);
}



int		myworkproc3(XtPointer client_data)
{
	if (level > (nlevels - 1))
		level = nlevels - 1;
	fill_pixmap();
	if (level == 0) {
		level = nlevels - 1;
		showpixmap();
		fprintf(stderr, "all done - ready for mouse input\n");
		return(1);
	} else {
		level--;
		return(0);
	}
}









