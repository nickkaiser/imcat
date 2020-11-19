/*
 * contour_stuff.c
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/*
#include <sys/types.h>
#include <sys/stat.h>
*/

#include "fits.h"

#define      MAX(x,y) (((x) > (y)) ? (x) : (y)) 
#define      MIN(x,y) (((x) < (y)) ? (x) : (y)) 

/*
#ifndef	S_IFIFO
#define S_IFIFO 0010000
#endif
*/

void	smcontour(float **f, int N1, int N2, int df, int nplus, int nminus, int comc, char *comv[])
{
	FILE	*fitsf, *smf;
	int	i, j;
	fitsheader	*fits;

	/* remove MAGIC pixels
	for (i = 0; i < N2; i++)
		for (j = 0; j < N1; j++)
			if (f[i][j] == MAGIC)
				f[i][j] = 0;

	/* write the macro file for sm */
	if (!(smf = fopen("xxxx.sm", "w")))
		error_exit("smcontour: failed to open xxxx.sm\n");
	fprintf(smf,"xxxx\n");
	fprintf(smf, "\tdev x11 -g 256x256\n");
	fprintf(smf, "\tredraw\n");
	fprintf(smf, "redraw\n");
	fprintf(smf, "\tdefine file_type fits\n");
	fprintf(smf, "\timage xxxx.fits\n");
	fprintf(smf, "\tlimits 0 %d 0 %d\n", N1 - 1, N2 - 1);
	fprintf(smf, "\tbox\n");
	fprintf(smf, "\txlabel j\n");
	fprintf(smf, "\tylabel i\n");
	if (nplus > 0) {
		fprintf(smf, "\tset lev = %d, %d, %d\n", 0, df * nplus, df);
		fprintf(smf, "\tlevels lev\n");
		fprintf(smf, "\tcontour\n");
	}
	if (nminus > 0) {
		fprintf(smf, "\tltype 2\n");
		fprintf(smf, "\tset lev = %d, %d, %d\n", -df * nminus, -df, df);
		fprintf(smf, "\tlevels lev\n");
		fprintf(smf, "\tcontour\n");
		fprintf(smf, "\tltype 0\n");
	}
	fclose(smf);

	/* write the xxxx.fits file */
	if (!(fitsf = fopen("xxxx.fits", "w")))
		error_exit("smcontour: failed to open xxxx.fits\n");
	fits = new2Dfitsheader(N1, N2, FLOAT_PIXTYPE);
	write2Dfloatimage(f, fits);
	fflush(fitsf);
	system("sm -m xxxx.sm < /dev/tty");
}




void	smprofile(float **f, int N1, int N2, int ic, int jc)
{
	FILE	*datf, *smf;
	float	*nsum, *fsum, *ffsum, *fav, *ffav, *sigma;
	int	bin, i, j, nbins, noplot = 0, imax, jmax;
	int	*ntot;

	/* we pass negative N1 to signify we just want to write profile data to stdout */
	if (N1 < 0) {
		noplot = 1;
		N1 = -N1;
	}

	/* figure out the most distant pixel */
	imax = MAX(ic, N2 - ic);
	jmax = MAX(jc, N1 - jc);
	/* set up the bins */
	nbins  = ceil(sqrt((double) (imax * imax + jmax * jmax)));
	ntot = (int *) calloc(nbins, sizeof(int));
	nsum = (float *) calloc(nbins, sizeof(float));
	fsum = (float *) calloc(nbins, sizeof(float));
	fav = (float *) calloc(nbins, sizeof(float));
	ffav = sigma = ffsum = (float *) calloc(nbins, sizeof(float));

	/* accumulate nsum, fsum. ffsum */
	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {
			bin = (int) floor(0.5 + sqrt((double) ((i - ic) * (i - ic) + (j - jc) * (j - jc))));
			ntot[bin]++;
			if (f[i][j] == FLOAT_MAGIC)
				continue;
			if (bin < 0 || bin >= nbins)
				continue;
			nsum[bin] += 1.0;
			fsum[bin] += (float) f[i][j];
			ffsum[bin] += (float) (f[i][j] * f[i][j]);
		}
	}
	for (bin = 0; bin < nbins; bin++) {
		if (nsum[bin]) {
			fav[bin] = fsum[bin] / nsum[bin];
			ffav[bin] = ffsum[bin] / nsum[bin];
			sigma[bin] = ffav[bin] - fav[bin] * fav[bin];
			sigma[bin] = (sigma[bin] > 0 ? sqrt(sigma[bin]) : 0.0);
		}
	}
	
	/* write the profile data to xxxx.dat or stdout*/
	if (noplot) {
		datf = stdout;
	} else {
		if (!(datf = fopen("xxxx.dat", "w")))
			error_exit("smprofile: failed to open xxxx.dat\n:");
	}
	fprintf(datf, "#      bin        <f>      sigma       fsum       npix\n");
	for (bin = 0; bin < nbins; bin++)
		fprintf(datf, "%10d %10.3e %10.3e %10.3e %10d\n", bin, fav[bin], sigma[bin], fsum[bin], ntot[bin]);
	fclose(datf);
	if (noplot)
		exit(0);

	/* write the sm macro file to xxxx.sm */
	if (!(smf = fopen("xxxx.sm", "w")))
		error_exit("smprofile: failed to open xxxx.sm\n:");
	fprintf(smf,"xxxx\n");
	fprintf(smf, "\tdev x11 -g 256x256\n");
	fprintf(smf, "\tredraw\n");
	fprintf(smf, "redraw\n");
	fprintf(smf, "\tdata xxxx.dat\n");
	fprintf(smf, "\tread {bin 1 fav 2 sigma 3}\n");
	fprintf(smf, "\tlimits bin fav\n");
	fprintf(smf, "\tbox\n");
	fprintf(smf, "\txlabel bin\n");
	fprintf(smf, "\tylabel <f>\n");
	fprintf(smf, "\tconnect bin fav\n");
	fprintf(smf, "\terrorbar bin fav sigma 2\n");
	fprintf(smf, "\terrorbar bin fav sigma 4\n");
	fclose(smf);

	/* fire up sm */
	system("sm -m xxxx.sm < /dev/tty");
	exit(0);
}



