#define	usage "\n\n\n\
NAME\n\
	makechart - overlay FITS image with object markers\n\
\n\
SYNOPSIS\n\
	makechart [option...] < catfile > fitsout\n\
		-r rname a	# draw or zap disk\n\
		-c 		# draw 16-32 pixel collar\n\
		-s sf		# scrunch factor (1)\n\
		-f image	# fitsfile \n\
		-m mask		# maskfile \n\
		-z		# zap circles or ellipses around object\n\
		-e ename efac	# draw an ellipse\n\
		-v circval	# paint circles etc this value (MAGIC)\n\
\n\
DESCRIPTION\n\
	\"makechart\" draws boxes or crosshairs round objects\n\
	Normally uses the image found in the cat header.\n\
	-f option forces it to use fitsfile\n\
	With -m option it draws rectangles from maskfile\n\
	which must be an 'lc' format catalogue containing\n\
	entries for a air of position vectors x1[2], x2[2]\n\
	for bottom-left and top-right corners respectively.\n\
	Use '-r rname a' to draw circle of radius a times the value of the object item\n\
	'rname', so use e.g. '-r rh 3' for circles 3 times the half-light radius.\n\
	Similarly, with -e option it will draw an ellipse with ellipticity\n\
	efac times the object ellipticity and with sqrt(a b) equal\n\
	to the radius as calculated above.\n\
\n\
	makechart uses the 'iostream' mechanism, so the 'fitsfile'\n\
	argument may be a fits generating command followed by '|'.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"		

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../utils/arrays.h"
#include "../../utils/stats_stuff.h"
#include "makechart.h"
#include "../../catlib/cat.h"
#include "../../utils/iostream.h"

#define MAGIC FLOAT_MAGIC

#define	CURS_MIN	8
#define CURS_MAX	16
#define	DEF_CURSOR_VAL	MAGIC	
#define	DEF_CIRCLE_VAL	MAGIC
#define	NW	0
#define	NE	1
#define	SE	2
#define	SW	3


static	float	cursval, circval;
	

main(int argc, char *argv[])	
{
	float 		**f;
	float		x, y, a;
	int		arg = 1, N1, N2, cursor = 1, circle = 0; 
	iostream	*fitsstream;
	FILE		*fitsf;
	fitsheader	*fits;
	cathead		*thecathead;
	object		*theobject;
	int		i, j, l, b, r, t, sf;
	double		ll, lb, lr, lt;
	char		maskfilename[64], *fitsfilename, maskline[1024], *rname, *ename;
	char		tempstring[128];
	FILE		*maskf;
	int		forcefitsfile = 0, pixtype;
	int		domask = 0, dozap = 0, docollar = 0, radiustype, doellipse, etype;
	float		efac, phi, e, emax, radius;
	double		Xf[2], R, RR, E[2];
	int		X[2];
	
	/* defaults */
	sf = 1;
	doellipse = 0;
	emax = 0.8;
	cursval = DEF_CURSOR_VAL;
	circval = DEF_CIRCLE_VAL;
	cursor = 1;
	circle = 0;
	
	while (arg < argc) {
		if (*argv[arg] != '-')
			error_exit(usage);
		switch (*(argv[arg++]+1)) {
			case 'r':
				cursor = 0;
				circle = 1;
				rname = argv[arg++];
				if (1 != sscanf(argv[arg++], "%f", &a))
					error_exit(usage);
				break;
			case 's':
				if (1 != sscanf(argv[arg++], "%d", &sf))
					error_exit(usage);
				break;
			case 'f':
				fitsfilename = argv[arg++];
				forcefitsfile = 1;
				break;
			case 'm':
				if (1 != sscanf(argv[arg++], "%s", maskfilename))
					error_exit(usage);
				domask = 1;
				break;
			case 'e':
				ename = argv[arg++];
				if (1 != sscanf(argv[arg++], "%f", &efac))
					error_exit(usage);
				doellipse = 1;
				break;
			case 'z':
				dozap = 1;
				break;
			case 'c':
				docollar = 1;
				break;
			case 'v':
				if (1 != sscanf(argv[arg++], "%f", &cursval))
					error_exit(usage);
				circval = cursval;
				break;
			default:
				error_exit(usage);
				break;
		}
	}


	thecathead = readcathead();
	theobject = newobject(thecathead);
	connectcatheadtoobject(theobject);
        setaddress(theobject, getobjectitemindex("x", theobject), Xf);
	if (doellipse) {
        	setaddress(theobject, getobjectitemindex(ename, theobject), E);
	}
	if (circle) {
        	setaddress(theobject, getobjectitemindex(rname, theobject), &R);
	}
	allocobjectcontents(theobject);
	

	if (forcefitsfile) {
		fitsstream = openiostream(fitsfilename, "r");
		if (!fitsstream) {
			error_exit("makechart : unable to open stream for input");
		}
		fitsf = fitsstream->f;
		/* fitsf = fopen(fitsfilename, "r"); */
	} else {
		fitsf = fopen(*((char **) getheaderitemaddress("fits_name", thecathead)), "r");
	}	
	if (!fitsf)
		error_exit("makechart: unable to open fits file for input");
	read2Dfloatimage(&f, &N1, &N2, &fits, fitsf);

	while (readobject(theobject)) {
		X[0] = (int) floor(0.5 + Xf[0] / sf);
		X[1] = (int) floor(0.5 + Xf[1] / sf);
		if (cursor) {
 			if (docollar) {
				drawcollar((int) X[1], (int) X[0], f, N1, N2);
			} else {
				drawcursor(Xf[0], Xf[1], f, N1, N2);
			}
		}
		if (circle) {
			if (R > 0) {
				radius = a * R;
			} else {
				radius = 0;
			}
			if (doellipse) {
				phi = 0.5 * atan2(E[1], -E[0]);
				e = efac * sqrt(E[0] * E[0] +  E[1] * E[1]);
				if (e > emax) {
					e = emax;
				}
				if (dozap) {
					zapellipse(Xf[0], Xf[1], radius, e, phi, f, N1, N2);
				} else {
					drawellipse(Xf[0], Xf[1], radius, e, phi, f, N1, N2);
				}
			} else {
				if (dozap) {
					zapcircle(Xf[0], Xf[1], radius, f, N1, N2);
				} else {
					drawcircle(Xf[0], Xf[1], radius, f, N1, N2);
				}
			}
		}
	}

	add_comment(argc, argv, fits);
	if (domask) {
		sprintf(tempstring, "lc -o x1 x2 < %s", maskfilename);
		if (!(maskf = popen(tempstring, "r")))
			error_exit("makechart: unable to open mask file\n");
		while (fgets(maskline, 1024, maskf)) {
			if (maskline[0] != '#') {
				sscanf(maskline, "%lf %lf %lf %lf", &ll, &lb, &lr, &lt);
				l = (int) floor(ll);
				b = (int) floor(lb);
				r = (int) floor(lr);
				t = (int) floor(lt);
				l = (l < 0 ? 0 : l);
				l = (l >= N1 ? N1 - 1 : l);
				r = (r < 0 ? 0 : r);
				r = (r >= N1 ? N1 - 1 : r);
				t = (t < 0 ? 0 : t);
				t = (t >= N2 ? N2 - 1 : t);
				b = (b < 0 ? 0 : b);
				b = (b >= N2 ? N2 - 1 : b);
/*				if (l < 0 || l >= N1 || r < 0 || r >= N1)
					continue;
				if (b < 0 || b >= N2 || t < 0 || t >= N2)
					continue;*/
				j = l;
				for (i = b; i <= t; i++)
					f[i][j] = cursval;
				j = r;
				for (i = b; i <= t; i++)
					f[i][j] = cursval;
				i = b;
				for (j = l; j <= r; j++)
					f[i][j] = cursval;
				i = t;
				for (j = l; j <= r; j++)
					f[i][j] = cursval;
			}
		}
	}
	write2Dfloatimage(f, fits);
	exit(0);
}




void	drawcursor(double x, double y, float **f, int N1, int N2)
{
	int d;
	int	i, j;

	i = (int) floor(y);
	j = (int) floor(x);

	
	for (d = CURS_MIN; d <= CURS_MAX; d++) {
		if (i + d >= 0 && i + d < N2 && j >= 0 && j < N1)
			f[i+d][j] = cursval;
		if (i - d >= 0 && i - d < N2 && j >= 0 && j < N1)
			f[i-d][j] = cursval;
		if (j + d >= 0 && j + d < N1 && i >= 0 && i < N2)
			f[i][j+d] = cursval;
		if (j - d >= 0 && j - d < N1 && i >= 0 && i < N2)
			f[i][j-d] = cursval;
	}
}

void	drawcollar(int io, int jo, float **f, int N1, int N2)
{
	int i, j, d, direction;
	
	drawcircle(io, jo, 16, f, N1, N2);
	drawcircle(io, jo, 32, f, N1, N2);
	for (d = floor(0.5 + 16 * sqrt(0.5)); d <= floor(0.5 + 32 * sqrt(0.5)); d++) {
		for (direction = 0; direction < 4; direction++) {
			switch (direction) {
				case NW:
					i = io + d;
					j = jo - d;
					break;
				case NE:
					i = io + d;
					j = jo + d;
					break;
				case SE:
					i = io - d;
					j = jo + d;
					break;
				case SW:
					i = io - d;
					j = jo - d;
					break;
			}
			if (i >= 0 && i < N2 && j >= 0 && j < N1)
				f[i][j] = cursval;
		}
	}
}


#define PI M_PI

void	drawcircle(double x, double y, float r, float **f, int N1, int N2)
{
	int	ii, jj, i, j;
	double phi, dphi;

	i = (int) floor(y);
	j = (int) floor(x);
	
	dphi = 1.0 / (5.0 * r);
	if (r < 1) {
		if (i >= 0 && i < N2 && j >= 0 && j < N1) {
			f[i][j] = circval;
		}
		return;
	}
	for (phi = 0; phi <= 2 * PI; phi += dphi) {
		ii = floor(y + r * cos(phi));
		jj = floor(x + r * sin(phi));
		if (ii >= 0 && ii < N2 && jj >= 0 && jj < N1)
			f[ii][jj] = circval;
	}
}


void	drawellipse(double x, double y, float r, float e, float phi0, float **f, int N1, int N2)
{
	int	ii, jj, i, j;
	double phi, dphi, a, b, I, J, c, s;
	
	i = (int) floor(y);
	j = (int) floor(x);
	
	dphi = 1.0 / (5.0 * r);
	if (r < 1) {
		if (i >= 0 && i < N2 && j >= 0 && j < N1) {
			f[i][j] = circval;
		}
		return;
	}
	a = r * (1 + e);
	b = r * (1 - e);
	c = cos(phi0);
	s = sin(phi0);
	for (phi = 0; phi <= 2 * PI; phi += dphi) {
		I = a * cos(phi);
		J = b * sin(phi);
		ii = floor(y + c * I - s * J);
		jj = floor(x + s * I + c * J);
		if (ii >= 0 && ii < N2 && jj >= 0 && jj < N1)
			f[ii][jj] = circval;
	}
}



void	zapcircle(double x, double y, float r, float **f, int N1, int N2)
{
	int	ii, jj, i, j, ir;
	double	rr, dx, dy;

	i = (int) floor(y);
	j = (int) floor(x);
	ir = (int) ceil(r);
	
	for (ii = i - ir; ii <= i + ir; ii++) {
		for (jj = j - ir; jj <= j + ir; jj++) {
			dx = jj + 0.5 - x;
			dy = ii + 0.5 - y;
			rr = dx * dx + dy * dy;
			if (rr > r * r)
				continue;
			if (ii >= 0 && ii < N2 && jj >= 0 && jj < N1)
				f[ii][jj] = MAGIC;
			}
	}
}

void	zapellipse(double x, double y, float r, float e, float phi0, float **f, int N1, int N2)
{
	int	ii, jj, i, j, ir;
	double 	phi, dphi, a, b, I, J, c, s, dx, dy, rx, ry, rr;
	
	i = (int) floor(y);
	j = (int) floor(x);
	
	dphi = 1.0 / (5.0 * r);
	if (r < 1) {
		if (i >= 0 && i < N2 && j >= 0 && j < N1) {
			f[i][j] = circval;
		}
		return;
	}
	a = r * (1 + e);
	b = r * (1 - e);
	c = cos(phi0);
	s = sin(phi0);
	ir = (int) ceil(a);
	for (ii = i - ir; ii <= i + ir; ii++) {
		for (jj = j - ir; jj <= j + ir; jj++) {
			dx = jj + 0.5 - x;
			dy = ii + 0.5 - y;
			rx = (c * dx - s * dy) / b;
			ry = (c * dy + s * dx) / a;
			rr = rx * rx + ry * ry;
			if (rr > 1) {
				continue;
			}
			if (ii >= 0 && ii < N2 && jj >= 0 && jj < N1) {
				f[ii][jj] = MAGIC;
			}
		}
	}
}


#undef PI

