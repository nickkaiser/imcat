#define	usage "\n\n\n\
NAME\n\
	minisky --- make a smooth local sky flat from cat of minima\n\
\n\
SYNOPSIS\n\
	minisky [option...] <catalogue >fits \n\
		-s n		# work with scrunch^n-ed image (2)\n\
		-r rf		# gaussian filter radius in real pixels (64)\n\
\n\
DESCRIPTION\n\
	\"minisky\" reads minima from a catalogue on stdin\n\
	It then creates two coarsely sampled images (chunky pixel side\n\
	= 2^n real pixels); f1 is number of minima falling in pixel\n\
	and f2 is sum of their fs values.  These are smoothed with\n\
	gaussian radius rf and then fsky = f2 / f1.\n\
	This will find low-frequency sky variations (as well as any\n\
	v low frequency and low surface brightness real objects.\n\
	Output image needs unscrunching to full resolution\n\
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
#include "../../imlib/filters.h"

#define MAGIC FLOAT_MAGIC

#define	SIGMA1		0.0
#define ALPHA		1.0
#define MAGICSUBSTITUTE	0.0
#define TINYVALUE	1.0e-6


main(int argc, char *argv[])	
{
	float 		**f1, **f2, scale = 1;
	int		arg = 1, Ng1, Ng2, nscrunch; 
	fitsheader	*fits;
	char		line[256], lcstring[64];
	int		ix, iy;
	float		rf, fs;
	FILE		*lcpipe;
	double		x, y;

	/* defaults */
	nscrunch = 2;
	rf = 64;
	strcpy(lcstring, "lc -P fits_size -o x fs");
	
	while (arg < argc) {
		if (*argv[arg] != '-')
			error_exit(usage);
		switch (*(argv[arg++]+1)) {
			case 'r':
				if (1 != sscanf(argv[arg++], "%f", &rf))
					error_exit(usage);
				break;
			case 's':
				if (1 != sscanf(argv[arg++], "%d", &nscrunch))
					error_exit(usage);
				break;
			default:
				error_exit(usage);
				break;
		}
	}
	
	if (!(lcpipe = popen(lcstring, "r")))
		error_exit("minisky: failed to open lc-pipe for input\n");
	fgets(line,  255, lcpipe);
	sscanf(line, "%d %d", &Ng1, &Ng2);
	while (nscrunch > 0) {
		Ng1 = ceil(0.5 * Ng1);
		Ng2 = ceil(0.5 * Ng2);
		scale *= 2;
		nscrunch--;
	}

	fits = new2Dfitsheader(Ng1, Ng2, FLOAT_PIXTYPE);
	add_comment(argc, argv, fits);
	
	/* make the two images */
	/* make them twice as big as necessary to avoid periodic bc's*/
	allocFloatArray(&f1, 2 * Ng1, 2 * Ng2); 
	allocFloatArray(&f2, 2 * Ng1, 2 * Ng2); 
	 
	while(fgets(line, 255, lcpipe)) {
		sscanf(line, "%lf %lf %f", &x, &y, &fs);
		ix = floor(x / scale);
		iy = floor(y / scale);
		if (ix >= 0 && ix < 2 * Ng1 && iy >= 0 && iy < 2 * Ng2 && fs != (float) MAGIC) {
			f1[iy][ix] += 1.0;
			f2[iy][ix] += fs;
		}
	}
	
	/* smooth the images */
	schecterfilter(f1, 2 * Ng1, 2 * Ng2, f1, SIGMA1, rf / scale, ALPHA, MAGICSUBSTITUTE);
	schecterfilter(f2, 2 * Ng1, 2 * Ng2, f2, SIGMA1, rf / scale, ALPHA, MAGICSUBSTITUTE);
 
	/* divide */
	for (iy = 0; iy < Ng2; iy++)
		for (ix = 0; ix < Ng1; ix++)
			f2[iy][ix] = (fabs(f1[iy][ix]) > TINYVALUE ? f2[iy][ix] / f1[iy][ix] : 0.0) ;
	write2Dfloatimage(f2, fits);
	exit(0);
}







