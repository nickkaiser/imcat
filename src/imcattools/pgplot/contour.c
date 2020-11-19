#define usage "\n\
NAME\n\
	contour --- create contour plot from fits image\n\
\n\
SYNOPSIS\n\
	contour [option....] < fitsin\n\
		-n nc		# number of contour intervals (10)\n\
		-f fmin fmax	# range of contour heights\n\
		-d device	# pgplot style device ('/xserve')\n\
		-t title	# text for title ('contour plot')\n\
		-l xlab ylab	# labels for axes ('x' 'y')\n\
		-w lwidth	# line width (1)\n\
		-H charheight	# character height (1)\n\
		-X X1 X2 Y1 Y2	# range for tick-labels\n\
		-j		# switch off justification\n\
		-a axis		# pgenv 'axis' value (0)\n\
		-g		# add gray-scale background image\n\
		-c colmap	# add color-mapped background image\n\
		-L x y text	# add a label\n\
		-W		# add a wedge for image range\n\
		-C		# output mouse-clicks\n\
\n\
DESCRIPTION\n\
	\"contour\" produces a contour plot from a fits image\n\
	using pgplot routines.\n\
\n\
	If fmin, fmax values are not specified these are\n\
	calculated from the input image.\n\
\n\
	It then draws (nc + 1) contours at levels\n\
	f - fmin + i * df, with df = (fmax - fmin) / nc.\n\
	By default, it produces output in an X-window on\n\
	the screen, but use -d option to specify alternative.\n\
\n\
	Use -c to display a colorised image - colmap can be 0,1,2.\n\
\n\
	Use -L option to add labels at arbitrary positions.\n\
\n\
	Use -W option to add a wedge showing the range of image values.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <float.h>

#include "cpgplot.h"

#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../utils/colormaps.h"
#include "labels.h"

#define MINFLOAT FLT_MIN
#define MAXFLOAT FLT_MAX

int		main(int argc, char *argv[])	
{
	int		arg = 1, x, y, N1, N2, nc, ic, lwidth;
	int		dograyscale, docolourimage;
	fitsheader	*fits;
	int		autorange;
	char		*device, defdevice[64] = "/xserve";
	float		**f, fmin, fmax, df, *ff, *flev;
	float		tr[6];
	char		*xlabel, *ylabel, *plotlabel, 
			defxlabel[128] = "x", defylabel[128] = "y", 
			defplotlabel[128] = "contour plot";
	float		x0, dx, y0, dy, x1, x2, y1, y2, X1, X2, Y1, Y2;
	int		just, axis;
	float		*cmapl, *cmapr, *cmapg, *cmapb, cmapcontra, cmapbright;
	int		c, cmapn, cmapindex, dowedge;
	label		*baselabel = NULL, *newlabel;
	float		charheight;
	int		usertickscale;
	int		getclicks;
	float		cursorx, cursory;
	char		cursevent[64];
	FILE		*opf;
	/* float		*xarray, *ylims;*/

	/* defaults */
	nc = 10;
	device = defdevice;
	autorange = 1;
	xlabel = defxlabel;
	ylabel = defylabel;
	plotlabel = defplotlabel;
	lwidth = 1;
	x0 = y0 = 0.0;
	dx = dy = 1.0;
	just = 1;
	axis = 0;
	dograyscale = docolourimage = 0;
	dowedge = 0;
	usertickscale = 0;
	charheight = 1;	
	getclicks = 0;

	while (arg < argc) {
                if (*argv[arg] != '-')
                        error_exit(usage);
                switch (*(argv[arg++]+1)) {
                        case 'n':
				if (1 != sscanf(argv[arg++], "%d", &nc))
					error_exit(usage);
                               break;
                        case 'f':
				if (1 != sscanf(argv[arg++], "%f", &fmin))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%f", &fmax))
					error_exit(usage);
				if (!(fmax > fmin))
					error_exit("contour: fmax must exceed fmin\n");
				autorange = 0;
                                break;
			case 'd':
				device = argv[arg++];
				break;
			case 't':
				plotlabel = argv[arg++];
				break;
			case 'l':
				xlabel = argv[arg++];
				ylabel = argv[arg++];
				break;
			case 'w':
				sscanf(argv[arg++], "%d", &lwidth);
				break;
			case 'H':
				sscanf(argv[arg++], "%f", &charheight);
				break;
			case 'x':
				error_exit("contour -x option disabled\n");
			case 'y':
				error_exit("contour -y option disabled\n");
			case 'X':
				sscanf(argv[arg++], "%f", &X1);
 				sscanf(argv[arg++], "%f", &X2);
				sscanf(argv[arg++], "%f", &Y1);
 				sscanf(argv[arg++], "%f", &Y2);
				usertickscale = 1;
				break;
			case 'j':
				just = 0;
				break;
			case 'a':
				sscanf(argv[arg++], "%d", &axis);
				break;
 			case 'g':
				dograyscale = 1;
				break;
 			case 'c':
				docolourimage = 1;
				sscanf(argv[arg++], "%d", &cmapindex);
				break;
			case 'L':
				/* create the new label */
				newlabel = (label *) calloc(1, sizeof(label));
				sscanf(argv[arg++], "%f", &(newlabel->x));
				sscanf(argv[arg++], "%f", &(newlabel->y));
				newlabel->text = argv[arg++];
				/* and install it at the front of the list */
				newlabel->next = baselabel;
				baselabel = newlabel;
				break;
			case 'W':
				dowedge = 1;
				break;
			case 'C':
				getclicks = 1;
				break;
                        default:
                                error_exit(usage);
                                break;
                }
        }
        
	
	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);
	ff = (float *) calloc(N1 * N2, sizeof(float));
	if (nc) {
		flev = (float *) calloc(nc + 1, sizeof(float));
	}

	if (autorange) {
		fmin = MAXFLOAT;
		fmax = -MAXFLOAT;
	}
	for (y = 0; y < N2; y++) {
		for (x = 0; x < N1; x++) {
			if (f[y][x] != FLOAT_MAGIC) {
				if (autorange) {
					if (f[y][x] > fmax) {
						fmax = f[y][x];
					}
					if (f[y][x] < fmin) {
						fmin = f[y][x];
					}
				}
				ff[N1 * y + x] = f[y][x];
			} else {
				ff[N1 * y + x] = 0.0;
			}
		}
	}
	if (autorange && (fmax == fmin)) {
		error_exit("contour: bad image: fmax = fmin\n");
	}
	if (nc) {
		df = (fmax - fmin) / nc;
		for (ic = 0; ic <= nc; ic++) {
			flev[ic] = fmin + ic * df;
		}
	}

	x1 = x0;
	x2 = x0 + N1 * dx;
	y1 = y0;
	y2 = y0 + N2 * dy;

	tr[0] = x0 - 0.5 * dx;
	tr[1] = dx;
	tr[2] = 0.0;
	tr[3] = y0 - 0.5 * dy;
	tr[4] = 0.0;
	tr[5] = dy;

	if (1 != cpgbeg(0, device, 1, 1))
		error_exit("contour: failed to open device\n");
	cpgslw(lwidth);
	cpgsch(charheight);
	cpgenv(x1, x2, y1, y2, just, -2);
	cpglab(xlabel, ylabel, plotlabel);
	cpgbbuf();

/*
	xarray = (float *) calloc(N1, sizeof(float));
	ylims = (float *) calloc(N1, sizeof(float));
	for (x = 0; x < N1; x++) {
		xarray[x] = x;
	}
	cpghi2d(ff, N1, N2, 1, N1, 1, N2, xarray, 1, 1.0, 0, ylims);
	cpgebuf();
	cpgend();
	exit(0);
*/

	if (dograyscale) {
		cpggray(ff, N1, N2, 1, N1, 1, N2, fmax, fmin, tr);
		if (dowedge) {
			cpgqch(&charheight);
			cpgsch((float) (0.4 * charheight));
			cpgwedg("RG", 1, 10, fmax, fmin, " ");
			cpgsch(charheight);
		}
	}
	if (docolourimage) {
		if (!getcolormap(&cmapl, &cmapr, &cmapg, &cmapb, &cmapcontra, &cmapbright, &cmapn, cmapindex)) {
			error_exit("contour: illegal colormap index\n");
		}
		cpgctab(cmapl, cmapr, cmapg, cmapb, cmapn, cmapcontra, cmapbright);
		cpgimag(ff, N1, N2, 1, N1, 1, N2, fmax, fmin, tr);
		if (dowedge) {
			cpgqch(&charheight);
			cpgsch((float) (0.4 * charheight));
			cpgwedg("RI", 1, 10, fmax, fmin, " ");
			cpgsch(charheight);
		}
	}
	if (nc) {
		cpgcont(ff, N1, N2, 1, N1, 1, N2, flev, nc, tr);
	}
	
        if (usertickscale) {
                cpgswin(X1, X2, Y1, Y2);
        }
	switch (axis) {
			case 0:
				cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);
				break;
			case 1:
				cpgbox("BCNSTA", 0.0, 0, "BCNSTA", 0.0, 0);
				break;
			case 2:
				cpgbox("BCNSTG", 0.0, 0, "BCNSTG", 0.0, 0);
				break;
			case -1:
				cpgbox("BC", 0.0, 0, "BC", 0.0, 0);
				break;
			case -2:
				cpgbox("", 0.0, 0, "", 0.0, 0);
				break;
               		case 10:
                        	cpgbox("BCLNST", 0.0, 0, "BCNST", 0.0, 0);
                        	break;
               		 case 20:
                        	cpgbox("BCNST", 0.0, 0, "BCLNST", 0.0, 0);
                        	break;
                	case 30:
                        	cpgbox("BCLNST", 0.0, 0, "BCLNST", 0.0, 0);
                        	break;
                	case 40:
                        	cpgtbox("BCNSTHZ", 0.0, 0, "BCNSTDZ", 0.0, 0);
                        	break;
                	case 42:
                        	cpgtbox("BCNSTHZG", 0.0, 0, "BCNSTDZG", 0.0, 0);
                        	break;
        	        case 43:
			        cpgtbox("BCNSTHZ", 0.0, 0, "BCNST", 0.0, 0);
	                	break;
	        	case 44:
	        	       cpgtbox("BCNST", 0.0, 0, "BCNSTHZ", 0.0, 0);
	        	       break;
  			default:
				error_exit("contour: illegal axis value\n");
	}
        if (usertickscale) {
                cpgswin(x1, x2, y1, y2);
        }

	/* now add the labels last */
	while (baselabel) {
		cpgtext(baselabel->x, baselabel->y, baselabel->text); 
		baselabel = baselabel->next;
	}
	cpgebuf();

	if (getclicks) {
		opf = popen("lc -C -N '1 2 x'", "w");
		while (1) {
			cpgcurs(&cursorx, &cursory, cursevent);
			switch (cursevent[0]) {
				case 'A':
					fprintf(opf, "%14.8g %14.8g\n", cursorx, cursory);
					break;
				case 'X':
					pclose(opf);
					cpgend();
					exit(0);
			}
		}
	}
	cpgend();

	exit(0);
}







