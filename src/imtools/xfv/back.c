/*
 * xfv - back.c
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
#include	"front.h"
#include	"back.h"
#include	"../../imlib/fits.h"
#include	"../../utils/error.h"
#include	"../../utils/colormaps.h"

#define		MAX_SHADES	4096

#define	MAX(x,y) (((x) > (y)) ? (x) : (y))
#define	MIN(x,y) (((x) < (y)) ? (x) : (y))

#define MAGIC SHORT_MAGIC

/* local globals */
static	int		ncolors;
static	int 		myscreen, depth;
static	GC		mygc;
static	Window		mywindow;
static	XImage		*theimage;

static Pixmap		pixmap;
static Colormap		cmap;
static Display 		*mydisplay;
static XColor		shade[MAX_SHADES];
static unsigned 	long	pixel[MAX_SHADES];
static unsigned int	n_shades;
static float		**f, gfmin, gfmax;
static int		N1, N2, blocksize, colormapindex, nframes; 
static fitsheader	*fits;
static float		timestep;


#define MAX_COLS 65535
#define	MAX(x,y) (((x) > (y)) ? (x) : (y))
#define	MIN(x,y) (((x) < (y)) ? (x) : (y))


/* readimheader  -  read the fits header, check dimensions and return size */
void	readimheader(int *n1, int *n2)
{
	fitscomment	*thecomment;

	fits = readfitsheader(stdin);
	if (fits->ndim != 3) { 
		error_exit("xfv: image dimensionality must be 3\n");
	}
	*n1 = N1 = fits->n[0];
	*n2 = N2 = fits->n[1];
	nframes = fits->n[2];
	if (thecomment = getcommentbyname("TIMESTEP", fits)) {
		timestep = getnumericvalue(thecomment);
	} else {
		timestep = 1.0 / nframes;
	}
	allocFloatArray(&f, N1, N2);
}


/*
 * set_shades(): (re)sets the colors in the part of the color map
 * we have reserved.  
 */
void	set_shades(void)
{
	int	i;
	long	longshade;
	float	*r, *g, *b;

	if (colormapindex >= 0) {
		getrgbfromcmap(&r, &g, &b, n_shades, colormapindex);
	}
	
	for (i = 0; i < n_shades; i++) {
		longshade = (long) (MAX_COLS * i / n_shades);
		shade[i].pixel = pixel[i];
		if (colormapindex < 0) {
			shade[i].red = shade[i].green = shade[i].blue = 
				(unsigned short) MAX(0, MIN(longshade, MAX_COLS));
		} else {
			shade[i].red   = (unsigned short) (MAX_COLS * MAX(0, MIN(r[i],1.0)));
			shade[i].green = (unsigned short) (MAX_COLS * MAX(0, MIN(g[i],1.0)));
			shade[i].blue  = (unsigned short) (MAX_COLS * MAX(0, MIN(b[i],1.0)));
		}
		shade[i].flags = DoRed | DoGreen | DoBlue;
	}
}

#undef MAX_COLS


void		alloc_shades(void)
{
	int		i, color, shadesgoal;
	XVisualInfo	info;
	Visual		*visual;
	int		class;
	unsigned long		plane_masks[1];

	mydisplay = XtDisplay(xfgetwidget(TOPLEVEL));
	myscreen = DefaultScreen(mydisplay);
	mywindow = DefaultRootWindow(mydisplay);	
	mygc = XCreateGC(mydisplay, mywindow, 0, 0);

	depth = DefaultDepth(mydisplay, myscreen);
	fprintf(stderr, "%d bit display\n", depth);
	shadesgoal = floor(0.5 + pow(2.0, (double) depth));
	/* fprintf(stderr, "shadesgoal = %d\n", shadesgoal); */
	cmap = DefaultColormap(mydisplay, myscreen);
	visual = DefaultVisual(mydisplay, myscreen);
	if (depth == 1) {
		fprintf(stderr, "this program won't work on a monochrome display\n");
		exit(0);
	}

	if (1) {
		class = PseudoColor;
		if (!XMatchVisualInfo(mydisplay, myscreen, depth, class, &info)) {
			fprintf(stderr, "no PseudoColor visual found\n");
			exit(0);
		}
		n_shades = shadesgoal;
		while (!XAllocColorCells(mydisplay, cmap, False, plane_masks, 0, pixel, n_shades)) {
			n_shades--;
			if (!n_shades) {
				error_exit("alloc_shades: unable to allocate any colors\n");
			}
		}
		fprintf(stderr, "allocated %d color cells\n", n_shades);
		set_shades();
		XStoreColors(mydisplay, cmap, shade, n_shades);
	} else {
		class = StaticColor;
		if (!XMatchVisualInfo(mydisplay, myscreen, depth, class, &info)) {
			fprintf(stderr, "no 8 bit StaticColor visual found\n");
			exit(0);
		}
		n_shades = 128;
		set_shades();
		for (i = 0; i < n_shades; i++) {
			if (!XAllocColor(mydisplay, cmap, &(shade[i]))) {
				fprintf(stderr, "alloc_shades: unable to allocate read only color\n");
			}
		}
	}
}



/*
 * creat_pixmap():  creates a pixmap
 */
void	create_pixmap(void) 
{
	pixmap = XCreatePixmap(mydisplay, mywindow, N1 * blocksize, N2 * blocksize, depth);
	fprintf(stderr, "memory reserved for %4d x %4d x %2d pixmap....\n", N1 * blocksize, N2 * blocksize, depth);
}

/*
 * fill_pixmap():  fills a pixmap 
 */
void	fill_pixmap(void) 
{
	int		x, y, ci, ix, iy;	
	unsigned long theshade;

	readfitsplane((void *) f, fits);	

	for (y = 0; y < N2; y++) {
		for (x = 0; x < N1; x++) {
			ci = MIN(n_shades - 1, MAX(0, (int) (n_shades * (f[y][x] - gfmin) / (gfmax - gfmin))));
			theshade = shade[ci].pixel;
			for (iy = 0; iy < blocksize; iy++) {
				for (ix = 0; ix < blocksize; ix++) {
					XPutPixel(theimage, x * blocksize + ix, y * blocksize + iy, theshade);
				}
			}
		}
	}
	XPutImage(mydisplay, pixmap, mygc, theimage, 0, 0, 0, 0, N1 * blocksize, N2 * blocksize);
}







/*
 * makes a XImage object which we will use to transfer
 * data to the pixmap
 */
void	makeXImage(void)
{
	unsigned int	format; 
	int		bitmap_pad;	
	int	width, height;
	
	width = blocksize * N1;
	height = blocksize * N2;

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





int		myworkproc1(XtPointer client_data)
{
	static int frame = 0;
	static	char string[64];

	if (frame >= fits->n[2]) {
		exit(0);
	} else {
		sprintf(string, "f=%5d\n", frame);
		setlabelstring(xfgetwidget(FLABEL), string);
		sprintf(string, "t=%12.3g\n", frame * timestep);
		setlabelstring(xfgetwidget(TLABEL), string);
		XFlush(mydisplay);
		fill_pixmap();
		showpixmap();
		frame++;
		return(0);
	}
}



void    showpixmap(void)
{
        Cardinal                n;
        Arg             arg[5];

        n = 0;
        XtSetArg(arg[n], XtNheight, N2 * blocksize);  n++;
        XtSetArg(arg[n], XtNwidth, N1 * blocksize); n++;
        XtSetArg(arg[n], XtNbitmap, pixmap); n++;
        XtSetValues(xfgetwidget(PICWIDGET), arg, n);
}



void	setblocksize(int b)
{
	blocksize = b;
}

void	setfrange(float ffmin, float ffmax)
{
	gfmin = ffmin;
	gfmax = ffmax;
}



void	setcolormapindex(int theindex)
{
	colormapindex = theindex;
}








