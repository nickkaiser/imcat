#define usage "\n\
NAME\n\
		listplot --- make 2-d scatter plot from a list of numbers\n\
\n\
SYNOPSIS\n\
		listplot [options....]\n\
			-d device	# pgplot style device ('/xserve')\n\
			-t title	# text for title ('xy-plot')\n\
			-l xlab ylab	# labels for axes ('x' 'y')\n\
			-w lwidth	# line width (1)\n\
			-x x1 x2	# left and right x-values\n\
			-y y1 y2	# bottom and top y-values\n\
			-j		# force justification\n\
			-a axis		# pgenv 'axis' value (0)\n\
			-m nmax		# max number of objects (100000)\n\
			-n symbol	# pgplot symbol (-1)\n\
\n\
DESCRIPTION\n\
		\"listplot\" reads a list of x,y pairs (one per line)\n\
		from stdin (ignorring lines beginning with '#') and\n\
		plots points.\n\
		Range of x, y values can be specified with -x, -y\n\
		options.  Otherwise they are calculated from input data.\n\
		By default, it produces output in an X-window on\n\
		the screen, but use -d option to specify alternative.\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>

#include "cpgplot.h"

#include "../../utils/error.h"

#define MARGIN 0.03

int		main(int argc, char *argv[])	
{
	int		arg = 1, lwidth, needxrange, needyrange;
	char		*device, defdevice[64] = "/xserve", line[1024];
	char		*xlabel, *ylabel, *plotlabel, 
			defxlabel[128] = "x", defylabel[128] = "y", 
			defplotlabel[128] = "xy-plot";
	float		xmin, xmax, ymin, ymax, x1, x2, y1, y2, *x, *y;
	int		just, axis, nmax, npts, nsym;


	/* defaults */
	device = defdevice;
	needxrange = 1;
	needyrange = 1;
	xlabel = defxlabel;
	ylabel = defylabel;
	plotlabel = defplotlabel;
	lwidth = 1;
	just = 0;
	axis = 0;
	nmax = 100000;
	nsym = -1;

	while (arg < argc) {
                if (*argv[arg] != '-')
                        error_exit(usage);
                switch (*(argv[arg++]+1)) {
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
			case 'x':
				sscanf(argv[arg++], "%f", &x1);
 				sscanf(argv[arg++], "%f", &x2);
				needxrange = 0;
				break;
			case 'y':
				sscanf(argv[arg++], "%f", &y1);
 				sscanf(argv[arg++], "%f", &y2);
				needyrange = 0;
				break;
			case 'j':
				just = 1;
				break;
			case 'a':
				sscanf(argv[arg++], "%d", &axis);
				break;
			case 'm':
				sscanf(argv[arg++], "%d", &nmax);
				break;
			case 'n':
				sscanf(argv[arg++], "%d", &nsym);
				break;
                        default:
                                error_exit(usage);
                                break;
                }
        }
        
	/* allocate space for data */
	x = (float *) calloc(nmax, sizeof(float));
	y = (float *) calloc(nmax, sizeof(float));

	/* read the data */
	npts = 0;
	while (fgets(line, 1024, stdin)) {
		if (line[0] == '#')
			continue;
		sscanf(line, "%f %f", x + npts, y + npts);
		if (!npts) {
			xmin = xmax = x[0];
			ymin = ymax = y[0];
		}
		if (needxrange) {
			if (x[npts] > xmax) xmax = x[npts];
			if (x[npts] < xmin) xmin = x[npts];
		}
		if (needyrange) {
			if (y[npts] > ymax) ymax = y[npts];
			if (y[npts] < ymin) ymin = y[npts];
		}
		npts++;
		if (npts == nmax)
			error_exit("listplot: too many points - increase nmax\n");
	}

	if (needxrange) {
		x1 = xmin - MARGIN * (xmax - xmin);
		x2 = xmax + MARGIN * (xmax - xmin);
	}
	if (needyrange) {
		y1 = ymin - MARGIN * (ymax - ymin);
		y2 = ymax + MARGIN * (ymax - ymin);
	}

	if (1 != cpgbeg(0, device, 1, 1))
		error_exit("contour: failed to open device\n");
	cpgslw(lwidth);
	cpgenv(x1, x2, y1, y2, just, axis);
	cpglab(xlabel, ylabel, plotlabel);
	cpgbbuf();
	cpgpt(npts, x, y, nsym);
	cpgebuf();
	cpgend();
	exit(0);
}







