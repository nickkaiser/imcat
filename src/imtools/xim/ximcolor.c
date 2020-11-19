/*
 * ximcolor.c
 *
 * color specific stuff for ximview X-image-viewer.
 * have separated this stuff so we chan chop and change color schemes
 *
 * this version does usual linear ramp for grayscale and in 3-color mode
 * in 2 color mode this uses "l-q" color model where q = ln(B/R) and
 * l = (R + B). 
 * 
 *	setcolorscheme(void)
 *	set_shades()
 *	color_index()
 */

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Intrinsic.h>
#include 	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<limits.h>
#include	"ximback.h"
#include	"ximcolor.h"
#include	"ximfront.h"
#include	"../../imlib/fits.h"
#include	"../../utils/error.h"

#define	MAX(x,y) (((x) > (y)) ? (x) : (y))
#define	MIN(x,y) (((x) < (y)) ? (x) : (y))
#define GREENFACTOR 2

/* for full caption */
/* #define CAPTION caption */
/* for postscript filename as caption */
#define CAPTION	psfilename
/* for no caption */
/* #define CAPTION	""*/

static	int	gfmin[MAX_COLORS], gfmax[MAX_COLORS];

/* externs */
extern	XColor		shade[];
extern	unsigned long	pixel[];
extern	unsigned int	n_shades;
extern	short		**f[MAX_FITS_DIM][MAX_LEVELS];
extern	Display 	*mydisplay;
extern	Colormap	cmap;
extern	int		N1[MAX_LEVELS], N2[MAX_LEVELS];
extern	char	caption[];

/* local globals */
static	int	ncc[MAX_COLORS], mul[MAX_COLORS], ngray, n_shades_used;
static	float	Fmin, Fmax, A, Y = GREENFACTOR, Q0, CFmax[3];	/* used for printing */
static	double		sliderval[MAX_COLORS];		
static	int	ncolors;


void		setcolorscheme(void)
{
	int	i = 1, sliderno;
	char	val[64];

	ncolors = xbgetvalue(NCOLORS);
	for (i = 0; i < ncolors; i++) {
		gfmin[i] = xfgetvalue(FMIN);
		gfmax[i] = xfgetvalue(FMAX);
	}
		
	switch (ncolors) {
		case	1:
			n_shades_used = ngray = n_shades;
			sprintf(val, "%.1f", newslidervalue(0));
			setslider(0, 1.0, "min", val);
			sprintf(val, "%.1f", newslidervalue(1));
			setslider(1, 0.0, "max", val);
			break;
		case	2:
			if (n_shades < 64)
				error_exit("setcolorscheme: too few shades\n");
			while (i * i <= n_shades)
				i++;
		    i--;
			n_shades_used = i * i;
			fprintf(stderr, "%d x %d colorbox\n", i, i);
			ncc[0] = i;
			ncc[1] = i;
			ncc[2] = 1;
			mul[0] = i;
			mul[1] = 1;
			mul[2] = 0;
			sprintf(val, "%.1f", gfmax[0]);
			setslider(0, 1.0, "fmax", val);
			setslider(1, 1.0 / A_MAX, "A", "1.0");
			setslider(2, 0.5, "q0", "0.0");
			break;
		case 	3:
			if (n_shades < 64)
				error_exit("setcolorscheme: too few shades\n");
			while (i * i * i <= n_shades)
				i++;
		    i--;
			n_shades_used = i * i * i;
			fprintf(stderr, "%d x %d x %d colorcube\n", i, i, i);
			ncc[0] = i;
			ncc[1] = i;
			ncc[2] = i;
			mul[0] = i * i;
			mul[1] = i;
			mul[2] = 1;
			sprintf(val, "%.1f", gfmax[0]);
			setslider(0, 1.0, "r", val);
			sprintf(val, "%.1f", gfmax[1]);
			setslider(1, 1.0, "g", val);
			sprintf(val, "%.1f", gfmax[2]);
			setslider(2, 1.0, "b", val);
			break;
		default:
			error_exit("setcolorscheme: bad ncolors\n");
			break;
	}

}

#define MAX_COLS 65535

/*
 * set_shades(): (re)sets the colors in the part of the color map
 * we have reserved.  Allows one to dynamically control the brighness
 * and contrast using the sliders.
 */
void	set_shades(void)
{
	int	ii, i, j, col[MAX_COLORS], ir, ig, ib, rgbmax;
	long	longshade, mycol[3];
	double	l, q, r, g, b, rgbsum;
	double	ramp, rampmin, rampmax;
	double	a, q0, *sval;

	sval = sliderval;
	switch (ncolors) {
		case 1:
			for (i = 0; i < n_shades; i++) {
				ramp = (double) i / (double) n_shades;
				rampmin = sval[0];
				rampmax = sval[1];
				if (rampmax != rampmin)
					longshade = (long) (MAX_COLS * (ramp - rampmin) / (rampmax - rampmin));
				else
					longshade = (long) (MAX_COLS * 0.5);
				shade[i].pixel = pixel[i];
				shade[i].red = shade[i].green = shade[i].blue = 
					(unsigned short) MAX(0, MIN(longshade, MAX_COLS));
				shade[i].flags = DoRed | DoGreen | DoBlue;
			}
			break;
		case 2:
			for (i = 0; i < ncc[0]; i++) {
				q = (double) (2 * i - ncc[0] + 1) / (double) ncc[0];
				for (j = 0; j < ncc[1]; j++) {
						l = (double) j / ((sval[0] + 0.01) * (ncc[1] - 1));
					a = A_MAX * sval[1];
					q0 = Q0_MAX * (2 * sval[2] - 1);
					r = 1.0;
					b = exp(a * (q - q0));
					g = GREENFACTOR * exp(0.5 * a * (q - q0));
					ii = i * mul[0] + j * mul[1];
					shade[ii].pixel = pixel[ii];
					rgbsum = r + g + b;
					if (rgbsum > 0) {
						ir = floor(0.5 + 3 * MAX_COLS * l * r / rgbsum);
						ig = floor(0.5 + 3 * MAX_COLS * l * g / rgbsum);
						ib = floor(0.5 + 3 * MAX_COLS * l * b / rgbsum);
					} else {
						ir = ig = ib = 0;
					}
					rgbmax = MAX(ir, MAX(ig, ib));
					if (rgbmax > MAX_COLS) {
						ir = floor(0.5 + ir * (double) MAX_COLS / (double) rgbmax);
						ig = floor(0.5 + ig * (double) MAX_COLS / (double) rgbmax);
						ib = floor(0.5 + ib * (double) MAX_COLS / (double) rgbmax);
					}
					shade[ii].red 	= (unsigned short) MIN(MAX_COLS, MAX(0, ir));
					shade[ii].green 	= (unsigned short) MIN(MAX_COLS, MAX(0, ig));
					shade[ii].blue 	= (unsigned short) MIN(MAX_COLS, MAX(0, ib));
					shade[ii].flags = DoRed | DoGreen | DoBlue;				
				}
			}
			break;
		case 3:
			for (col[0] = 0; col[0] < ncc[0]; col[0]++) {
				mycol[0] = col[0] * MAX_COLS / ((ncc[0] - 1) * sval[0]);
				for (col[1] = 0; col[1] < ncc[1]; col[1]++) {
					mycol[1] = col[1] * MAX_COLS / ((ncc[1] - 1) * sval[1]);
					for (col[2] = 0; col[2] < ncc[2]; col[2]++) {
						mycol[2] = col[2] * MAX_COLS / ((ncc[2] - 1) * sval[2]);
						i = col[2] * mul[2] + col[1] * mul[1] + col[0] * mul[0];
						shade[i].pixel = pixel[i];
						shade[i].red 	= (unsigned short) MIN(MAX_COLS, mycol[0]);
						shade[i].green 	= (unsigned short) MIN(MAX_COLS, mycol[1]);
						shade[i].blue 	= (unsigned short) MIN(MAX_COLS, mycol[2]);
						shade[i].flags = DoRed | DoGreen | DoBlue;				
					}
				}
			}
			break;
		default:
			error_exit("set_shades: bad ncolors\n");
	}
	XStoreColors(mydisplay, cmap, shade, n_shades_used);
}

#undef MAX_COLS



int	color_index(int level, int i, int j)
{
	int 	color, theindex, ii[MAX_COLORS], igray, df[MAX_COLORS], iq, il;
	double	l, q;

	theindex = 0;
	switch (ncolors) {
		case 1:
			df[0] = gfmax[0] - gfmin[0];
			igray = (f[0][level][i][j] - gfmin[0]) * ngray / df[0];
			theindex = MAX(0, MIN(igray, ngray - 1));
			break;
		case 2:
				l = (double) (f[0][level][i][j] + f[1][level][i][j]);
				if (f[0][level][i][j] && f[1][level][i][j])
					q = log(fabs((double) f[1][level][i][j] / (double) f[0][level][i][j]));
				else
					q = 0;
				il = floor(0.5 + ncc[0] * l / gfmax[0]);
				il = MAX(0, MIN((ncc[0] - 1), il));
				iq = ncc[1] * (1 + q) / 2;
				iq = MAX(0, MIN((ncc[1] - 1), iq));
				theindex = mul[0] * iq + mul[1] * il;
			break;
		case 3:
			for (color = 0; color < 3; color++) {
				df[color] = gfmax[color] - gfmin[color];
				ii[color] = (f[color][level][i][j] - gfmin[color]) * ncc[color] / df[color];
				theindex += mul[color] * MAX(0, MIN(ii[color], ncc[color] - 1));
			}
			break;
		default:
			break;
	}
	return (theindex);
}



float	newslidervalue(int index)
{
	switch (ncolors) {
		case 1:
			switch (index) {
				case 0:
					Fmin = gfmin[0] + sliderval[index] * (gfmax[0] - gfmin[0]);
					break;
				case 1:
					Fmax = gfmin[0] + sliderval[index] * (gfmax[0] - gfmin[0]);
					break;
				case 2:
					break;
				default:
					error_exit("newslidervalue: bad index\n");
					break;
			}
			return ((float) (gfmin[0] + sliderval[index] * (gfmax[0] - gfmin[0])));
			break;
		case 2:
			switch (index) {
				case 0:
					Fmax = sliderval[index] * gfmax[index];
					return (Fmax);
					break;
				case 1:
					A = sliderval[index] * A_MAX;
					return (A);
					break;
				case 2:
					Q0 = Q0_MAX * (sliderval[index] - 0.5);
					return (Q0);
					break;
				default:
					error_exit("newslidervalue: bad ncolors\n");
					break;
			}
			break;
		case 3:
			CFmax[index] = sliderval[index] * gfmax[index];
			return (sliderval[index] * gfmax[index]);
			break;
		default:
			error_exit("newslidervalue: bad ncolors\n");
	}
}

/* printing stuff.... */

static int clevel, io, jo;

void	colorsetlevel(int level)
{
	clevel = level;
}

void	printall(int level, char *psfilename)
{
	FILE	*opf;
	
	io = jo = 0;
	fprintf(stderr, "printing full view at level %d to %s\n", level, psfilename);
	colorsetlevel(level);
	opf = fopen(psfilename, "w");
	if (!opf)
		error_exit("printall: failed to open ps file for output\n");
	set_print_opf(opf);
	if (ncolors == 1)
		print_im(N1[level], N2[level], 1, gray, CAPTION);
	else
		print_im(N1[level], N2[level], 3, gray, CAPTION);
	fclose(opf);
}


void	printselection(int level, int i1, int j1, int i2, int j2, char *psfilename)
{
	FILE	*opf;
	int	xN1, xN2, mag = 1, lev = level;
	
	while (lev > 0) {
		mag *= 2;
		lev --;
	}
	xN2 = (i2 - i1) / mag;
	xN1 = (j2 - j1) / mag;
	io = i1 / mag;
	jo = j1 / mag;
	fprintf(stderr, "printing %d x %d selection at level %d to %s\n", xN1, xN2, level, psfilename);
	colorsetlevel(level);
	opf = fopen(psfilename, "w");
	if (!opf)
		error_exit("printselection: failed to open ps file for output\n");
	set_print_opf(opf);
	if (ncolors == 1)
		print_im(xN1, xN2, 1, gray, CAPTION);
	else
		print_im(xN1, xN2, 3, gray, CAPTION);
	fclose(opf);
}


double	gray(int i, int j, int color)
{
	static double	shade[3], shadesum, maxshade, q, l;
	static int	lasti, lastj, col;

	i += io;
	j += jo;
	if ((ncolors == 1 && color > 0) || (color > 2))
		error_exit("gray: bad color\n");
	if (i < 0 || i >= N2[clevel] || j < 0 || j >= N1[clevel])
		return(0.0);
	switch (ncolors) {
		case 1:
			return (1 - (double) (f[0][clevel][i][j] - Fmin) / (double) (Fmax - Fmin));
			break;
		case 2:
			if (i == lasti && j == lastj)		/* we already have it */
				return(1 - shade[color]);
			lasti = i;
			lastj = j;
			if (f[0][clevel][i][j] > 0 && f[1][clevel][i][j] > 0)
				q = log((double) f[1][clevel][i][j] / (double) f[0][clevel][i][j]);
			else
				q = 0.0;
			shade[0] = 1.0;
			shade[1] = GREENFACTOR * exp(0.5 * A * (q - Q0));
			shade[2] = exp(A * (q - Q0));
			shadesum = shade[0] + shade[1] + shade[2];
			l = f[0][clevel][i][j] + f[1][clevel][i][j];
			for (col = 0; col < 3; col++)
				shade[col] *= 3 * l / (shadesum * Fmax);
			maxshade = MAX(shade[0], MAX(shade[1], shade[2]));
			if (maxshade > 1.0)
				for (col = 0; col < 3; col++)
					shade[col] /= maxshade;
			return(1 - shade[color]);
			break;
		case 3:
			return (1 - (double) f[color][clevel][i][j] / (double) CFmax[color]);
			break;
		default:
			error_exit("gray: bad ncolors\n");
			break;
	}
}


void	setsliderval(int index, double theval)
{
	sliderval[index] = theval;
}






