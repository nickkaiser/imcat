/*
 * ximcolor.c
 *
 * color specific stuff for ximview X-image-viewer.
 * have separated this stuff so we chan chop and change color schemes
 *
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
#include	"back.h"
#include	"color.h"
#include	"front.h"
#include	"../../imlib/fits.h"
#include	"../../utils/error.h"

#define	MAX(x,y) (((x) > (y)) ? (x) : (y))
#define	MIN(x,y) (((x) < (y)) ? (x) : (y))

/* for full caption */
/* #define CAPTION caption */
/* for postscript filename as caption */
#define CAPTION	psfilename
/* for no caption */
/* #define CAPTION	""*/

static	int	fmin, fmax;

/* externs */
extern	XColor		shade[];
extern	unsigned long	pixel[];
extern	unsigned int	n_shades;
extern	Display 	*mydisplay;
extern	Colormap	cmap;
extern	int		N1, N2;
extern float		**f, fmin, fmax;
extern	char		caption[];

/* local globals */
static	int	ncc[MAX_COLORS], mul[MAX_COLORS], ngray, n_shades_used;



#define MAX_COLS 65535

/*
 * set_shades(): (re)sets the colors in the part of the color map
 * we have reserved.  
 */
void	set_shades(void)
{
	int	i;
	long	longshade, mycol[3];
	double	ramp, rampmin, rampmax;

	for (i = 0; i < n_shades; i++) {
		ramp = (double) i / (double) n_shades;
		rampmin = 0;
		rampmax = 255;
		if (rampmax != rampmin)
			longshade = (long) (MAX_COLS * (ramp - rampmin) / (rampmax - rampmin));
		else
			longshade = (long) (MAX_COLS * 0.5);
		shade[i].pixel = pixel[i];
		shade[i].red = shade[i].green = shade[i].blue = 
			(unsigned short) MAX(0, MIN(longshade, MAX_COLS));
		shade[i].flags = DoRed | DoGreen | DoBlue;
	}
	XStoreColors(mydisplay, cmap, shade, n_shades_used);
}

#undef MAX_COLS








