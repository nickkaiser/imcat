#define usage "\n\n\n\
NAME\n\
	mask - set abnormal pixels in sky-flat to MAGIC\n\
\n\
SYNOPSIS\n\
	mask [options] < fits_in > fits_out \n\
		-b inner outer 	# box for local median (3 5)\n\
		-f dfcrit	# critical (f - fmedian) value (100)\n\
\n\
DESCRIPTION\n\
	\"mask\" reads fits_in (typically a sky_flat) \n\
	calculates a local median and sets to MAGIC\n\
	any pixels with abs(f - f_local_median) > dfcrit \n\
	Finally any neighbours of MAGIC pixels are also set to MAGIC.\n\
	The MAGIC property will be inherited when the result is\n\
	used to sky flatten an image\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/arrays.h"

#define MAGIC FLOAT_MAGIC

int		floatcmp();

int		main(int argc, char *argv[])	
{
	int	arg = 1;
	int	i, j, di, dj, ii, jj, N1, N2, inner, outer;
	int	control_size, control_pts;
	fitsheader	*fits;
	float	**f, *fcontrol, **fdiff, dfcrit;
	double	fmedian;

	/* defaults */
	dfcrit = 100;
	inner = 3;
	outer = 5;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch (argv[arg++][1]) {
			case 'b':
				if (1 != sscanf(argv[arg++], "%d", &inner))
					error_exit(usage);
				if (1 != sscanf(argv[arg++], "%d", &outer))
					error_exit(usage);
				break;
			case 'f':
				if (1 != sscanf(argv[arg++], "%f", &dfcrit))
					error_exit(usage);
				break;
			default:
				error_exit(usage);
		}
	}
	
	control_size = (2 * outer + 1) * (2 * outer + 1) -
		(2 * inner - 1) * (2 * inner - 1);
	fcontrol = (float *) calloc(control_size, sizeof(float));

	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);
	allocFloatArray(&fdiff, N1, N2);
	

	for (i = 0; i < N2; i++) {
	  for (j = 0; j < N1; j++) {
	    control_pts = 0;
	    for (di = -outer; di <= outer; di++) {
	      if (di > -inner && di < inner)
		continue;
	      ii = i + di;
	      if (ii < 0 || ii >= N2)
		continue;
	      for (dj = -outer; dj <= outer; dj++) {
	        if (dj > -inner && dj < inner)
		  continue;
		jj = j + dj;
		if (jj < 0 || jj >= N1)
		  continue;
		fcontrol[control_pts++] = f[ii][jj];
 	      }
            }
	    qsort(fcontrol, control_pts, sizeof(float), floatcmp);
	    fmedian = 0.5 * (fcontrol[control_pts / 2] +
				fcontrol[control_pts / 2 - 1]);
	    fdiff[i][j] = fabs(f[i][j] - fmedian);
	  }
 	}

	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {
 			if (fdiff[i][j] > dfcrit) {
				f[i][j] = MAGIC;
			}
			fdiff[i][j] = f[i][j];
		}
	}

	for (i = 1; i < N2 - 1; i++) {
		for (j = 1; j < N1 - 1; j++) {
 			if (fdiff[i][j] == MAGIC) { 
				f[i-1][j-1] = f[i-1][j] = f[i-1][j+1] = MAGIC;
				f[i+1][j-1] = f[i+1][j] = f[i+1][j+1] = MAGIC;
				f[i][j-1] = f[i][j+1] = MAGIC;
			}
		}
	}

	add_comment(argc, argv, fits);
	write2Dfloatimage(f, fits);
	exit(0);	
}



int		floatcmp(float *s1, float *s2) 
{
	if (*s1 < *s2)
		return (-1);
	else if (*s1 > *s2)
		return (1);
	else
		return (0);
}

