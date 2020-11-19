/*
 * iisdisplay.c
 *
 */

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "iis.h"
#include "iisdisplay.h"

static	short		*data = NULL;
static	FILE		*iispipe, *syspipe;

#define	FMAX			200
#define	IIS_GREEN		201
#define IIS_BLACK		202
#define	IIS_WHITE		203
#define IIS_RED			204
#define IIS_BLUE		206
#define IIS_YELLOW		207

int	openiispipe(void)
{
	char	pipename[256];

	syspipe = popen("echo $HOME/dev/saopipe", "r");
	fscanf(syspipe, "%s", pipename);
	pclose(syspipe);
/*
	printf("pipename = %s\n", pipename);
	exit(0);
*/	
	if (!(iispipe = fopen(pipename, "w"))) {
		fprintf(stderr, "openiispipe: cannot open pipe\n");
		exit(-1);
	}
	return(1);
}

int	iisdisplay(float **f, int N1, int N2, float fmin, float fmax)
{
	unsigned	short	hdr[8];
	short		dataval;
	int		ix, iy, x, y, xsize, ysize, indx, checksum;
	float		ff, fscale = 1, scale, xscale, yscale;

	/* size of the image as it will apppear in the saoimage window */
	xsize = 512;
	ysize = 512;

	/* figure out the spatial scaling to get image into 512 x 512 box */
	xscale = N1 / (float) xsize;
	yscale = N2 / (float) ysize;
	scale = (xscale > yscale ? xscale : yscale);


	/* get max, min */
	if (fmin == 0.0 && fmax == 0.0) {
		fmin = (float) SHRT_MAX;
		fmax = (float) SHRT_MIN;
		for (y = 0; y < N2; y++) {
			for (x = 0; x < N1; x++) {
				ff = f[y][x];
				if ((ff != SHRT_MIN) && (ff != SHRT_MAX)) {
					fmin = (ff < fmin ? ff : fmin);
					fmax = (ff > fmax ? ff : fmax);
				}
			}
		}
	}
	if (fmax != fmin)
		fscale = FMAX / (fmax - fmin);


	if (fmax != fmin)
		fscale = FMAX / (fmax - fmin);
	else
		fscale = 1.0;

	/* create data */
	if (!data) {
		data = (short *) calloc(xsize, sizeof(short));
	}

	for (y = 0; y < ysize; y++) {
		/* create header */
		hdr[TRANSFER_ID] = IWRITE;
		hdr[THING_COUNT] = -xsize;
		hdr[SUB_UNIT] = REFRESH;
		hdr[CHECK_SUM] = 0;
		hdr[X_REGISTER] = ADVXONTC;
		hdr[Y_REGISTER] = ADVYONXOV + y;
		hdr[Z_REGISTER] = CHAN2;
		hdr[T_REGISTER] = ALLBITPL;
		checksum = 0;
		for (indx = 0; indx < 8; indx++) {
			checksum += hdr[indx];
		}
		hdr[CHECK_SUM] = 0177777 - (unsigned short) checksum;

		fwrite(hdr, sizeof(short), 8, iispipe);
		iy = (int) (scale * (ysize - 1 - y));
		for (x = 0; x < xsize; x++) {
			ix = (int) (scale * x);
			if (ix < 0 || ix >= N1 || iy < 0 || iy >= N2) {
				data[x] = 0;
			} else {
				if (f[iy][ix] == FLOAT_MAGIC) {
						data[x] = IIS_YELLOW;
				} else {
						data[x] = (int) (fscale * (f[iy][ix] - fmin));
						data[x] = IIS_YELLOW;
				}
			}
		}
		fwrite(data, sizeof(short), xsize, iispipe);
	}
	return (1);
}
