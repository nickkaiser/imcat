/*
 * register.c
 *
 * determine rotation, scale and offset for pair of lists of coords
 */


#include <stdio.h>
#include <math.h>
#include "../utils/arrays.h"
#include "../fftlib/myfft.h"
#include "../utils/error.h"
#include "../imlib/filters.h"

#define usage "usage: register xylist x'y'list [options...]\n\
	read 2 lists of coords: (x,y), (x',y') and determine\n\
	scale, rotation and translation that maps (x',y') => (x,y)\n\
		x = a (x' cos - y' sin) + x0\n\
		y = a (x' sin + y' cos) + y0\n\
	It works by cross-correlating images of wrapped lnr, phi\n\
	values for pairs of objects and locating peak to determine\n\
	a, phi.  We then rotate and scale x' y', generate images\n\
	of x,y positions and cross correlate these to get x0, y0.\n\
	Outputs x0, y0, a, phi.\n\
	Options are:\n\
		-i imsize	# internal image size (512)\n\
		-p phi1 phi2	# range for phi wrapping (0, PI)\n\
		-a a1 a2	# find peaks only in this range of scale\n\
		-d dlnr		# range of lnr for wrapping (1.0)\n\
\n\n"

#define	PI M_PI
#define MAX_OBJECTS 100000
#define BIG_NEG		-1.e10
#define BIG_POS		1.e10;

float   mexicanfilterfunction(float ki, float kj);

main(int argc, char *argv[])
{
	char	line[1024];
	double	rr, dx, dy, phi, dphi, phi1, phi2, lnr, dlnr, a1, a2, temp;
	int	arg, i, ix, iy, imsize, nobj[2], i0, i1; 
	FILE	*ipf[2], *pipe;
	float	**n[2], **theccf, sigma1, sigma2;
	double	ccfmax, ccfabsmax, phimax, amax, X, Y, c, s, x1, y1, x2, y2;
	double	x[2][MAX_OBJECTS], y[2][MAX_OBJECTS];
	fft_type	fk[2];

	/* defaults */
	phi1 = 0.0;
	phi2 = PI;
	dlnr = 1.0;
	imsize = 512;
	a1 = 0.1;
	a2 = 10.0;
	sigma1 = 1;
	sigma2 = 2;

	if (argc < 3) {
		fprintf(stderr, "%s", usage);
		exit(-1);
	}

	/* open the list files */
	if (!(ipf[0] = fopen(argv[1],"r")) || !(ipf[1] = fopen(argv[2],"r")))
		error_exit("register: failed to open coordinate list files\n");

	/* parse args */
	arg = 3;
	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch (argv[arg++][1]) {
			case 'i':
				sscanf(argv[arg++], "%d", &imsize);
				break;
			case 'p':
				sscanf(argv[arg++], "%lf", &phi1);
				sscanf(argv[arg++], "%lf", &phi2);
				break;
			case 'a':
				sscanf(argv[arg++], "%lf", &a1);
				sscanf(argv[arg++], "%lf", &a2);
				if (a2 < a1) {
					temp = a2;
					a2 = a1;
					a1 = temp;
				}
				break;
			case 'd':
				sscanf(argv[arg++], "%lf", &dlnr);
				break;
			default:
				error_exit(usage);
				break;
		}
	}
	dphi = phi2 - phi1;

	/* allocate image arrays, fft's */
	allocFloatArray(&theccf, imsize, imsize);
	for (i = 0; i < 2; i++) {
		allocFloatArray(&(n[i]), imsize, imsize);
		alloc_fft(&(fk[i]), imsize, imsize);
	}

  
	/* read the objects */
	for (i = 0; i < 2; i++) {
		nobj[i] = 0;
		while (fgets(line, 1024, ipf[i])) {
			if (line[0] == '#' || line[0] == '\0')
				continue;
			sscanf(line, "%lf %lf", &(x[i][nobj[i]]), &y[i][nobj[i]]);
			if (nobj[i]++ >= MAX_OBJECTS)
				error_exit("register: too many objects\n");
		}
	}

	/* make images of pairs */
	for (i = 0; i < 2; i++) {
		for (i0 = 0; i0 < nobj[i]; i0++) {
			for (i1 = 0; i1 < nobj[i]; i1++) {
				dx = x[i][i1] - x[i][i0];
				dy = y[i][i1] - y[i][i0];
				rr = dx * dx + dy * dy;
				if (rr == 0.0)
					continue;
				lnr = 0.5 * log(rr);
				phi = atan2(dy, dx) - phi1;
				ix = imsize * (phi/dphi - floor(phi/dphi));
				iy = imsize * (lnr/dlnr - floor(lnr/dlnr));
				n[i][iy][ix] += 1.0;
			}
		}
	}

	/* do the crosscorrelation */
	forward_fft(n[0], imsize, imsize, fk[0]);
	forward_fft(n[1], imsize, imsize, fk[1]);
	ccf(fk[0], fk[1], imsize, imsize, theccf, imsize / 2, imsize / 2);
	mexicanfilter(theccf, imsize, imsize, theccf, sigma1, sigma2, 0);

	/* find the cross correlation peak */
	ccfmax = ccfabsmax = 0.0;
	for (iy = 0; iy < imsize; iy++) {
		for (ix = 0; ix < imsize; ix++) {
			if (theccf[iy][ix] > ccfabsmax) {
				ccfabsmax = theccf[iy][ix];
			} 
			if (theccf[iy][ix] > ccfmax) {
				amax = exp(dlnr * (iy - imsize / 2) / imsize);
				phimax = phi1 + dphi * (ix - imsize / 2) / imsize; 
				if (amax > a1 && amax < a2) {
					ccfmax = theccf[iy][ix];
				}
			} 
		}
	}
	if (ccfabsmax != ccfmax)
		error_exit("register: found higher peak outside of a1-a2\n");
	phimax *= -1;
	amax = 1.0 / amax;

	/* now we scale and rotate the 2nd list of coords */
	c = cos(phimax);
	s = sin(phimax);
	for (i1 = 0; i1 < nobj[1]; i1++) {
		X = amax * (c * x[1][i1] - s * y[1][i1]);
		Y = amax * (s * x[1][i1] + c * y[1][i1]);
		x[1][i1] = X;
		y[1][i1] = Y;
	}
	
	/* now we figure out the bounding box */
	x1 = y1 = BIG_POS;
	x2 = y2 = BIG_NEG;
	for (i = 0; i < 2; i++) {
		for (i0 = 0; i0 < nobj[i]; i0++) {
			if (x[i][i0] > x2)
				x2 = x[i][i0];
			if (x[i][i0] < x1)
				x1 = x[i][i0];
			if (y[i][i0] > y2)
				y2 = y[i][i0];
			if (y[i][i0] < y1)
				y1 = y[i][i0];
		}
	}

	/* now generate the two images of x,y coords */
	for (iy = 0; iy < imsize; iy++) {
		for (ix = 0; ix < imsize; ix++) {
			n[0][iy][ix] = n[1][iy][ix] = 0.0;
		}
	}
	for (i = 0; i < 2; i++) {
		for (i0 = 0; i0 < nobj[i]; i0++) {
			ix = floor(imsize * (x[i][i0] - x1) / (x2 - x1));
			iy = floor(imsize * (y[i][i0] - y1) / (y2 - y1));
			if (ix >= 0 && ix < imsize && iy >= 0 && iy < imsize)
				n[i][iy][ix] += 1.0;
		}
	}

	/* now we cross correlate and smooth */
	forward_fft(n[0], imsize, imsize, fk[0]);
	forward_fft(n[1], imsize, imsize, fk[1]);
	ccf(fk[0], fk[1], imsize, imsize, theccf, imsize / 2, imsize / 2);
	mexicanfilter(theccf, imsize, imsize, theccf, sigma1, sigma2, 0);

/*
	set_output_pixtype(-32);
	set_fits_opf(fopen("temp1.fits", "w"));
	fwrite_fits(n[0], imsize, imsize, 0, argv);
	set_fits_opf(fopen("temp2.fits", "w"));
	fwrite_fits(n[1], imsize, imsize, 0, argv);
	set_fits_opf(fopen("temp3.fits", "w"));
	fwrite_fits(theccf, imsize, imsize, 0, argv);

*/

	/* and find the peak */
	ccfmax = 0.0;
	for (iy = 0; iy < imsize; iy++) {
		for (ix = 0; ix < imsize; ix++) {
			if (theccf[iy][ix] > ccfmax) {
				dx = -(x2 - x1) * (ix - imsize / 2) / imsize; 
				dy = -(y2 - y1) * (iy - imsize / 2) / imsize; 
				ccfmax = theccf[iy][ix];
			} 
		}
	}

	/* and we're all done */
	fprintf(stdout, "#       dx         dy          a        phi\n");
	fprintf(stdout, "%10.3lf %10.3lf %10.6lf %10.6lf\n", dx, dy, amax, phimax);	
}


 