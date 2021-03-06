/**
 **
 ** routines for calulating statistics of fields:
 **
 ** mostly straightforward 'cept for the mode & sigma where we first
 ** make a crude estimate of the mode as the most frequently
 ** occuring value and estimate sigma from sub-mode f values. We then
 ** refine this by making a least squares fit to the mode � sigma / 2 region
 **
 ** aug 4th '92:
 ** revised fdo_stats
 **	now estimates quartiles and median from a random sample
 ** 	calculates sigma from width of inner 25 percent range
 ** Fri Aug  7 09:29:37 EDT 1992
 **		do_stats also converted to sampling technique
 **		mode determined by smoothed histogram method
 **		incorporates findmode() - formerly in analyse.c
 **		general purpose routine liststats() to return
 **		mode, medians, quartiles from an unsorted list.
 **
 ** Tue Nov 25 11:21:48 HST 2003 fixed the mean bug in fdo_stats
 **/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
 
#include	"stats_stuff.h"
#include	"error.h"
#include	"ran1.h"
#include	"../imlib/fits.h"

#define	good(x)		(x != SHRT_MIN && x != SHRT_MAX && x != SHORT_MAGIC)
#define	fgood(x)	(x != FLOAT_MAGIC)

double	drand48();

#define DEFSAMPLEFRAC 0.1
#define HUGE_VAL	1.0e20
#define MAX_SAMPLE_SIZE	10000





void		fdo_stats(float **f, int N1, int N2, int margin, fstatsrec *srec)
{
	int		i, j, isample, ngood, badpix, step;
	int		actualsize, samplesize, gotagoodone;
	float		fmin, fmax, samplefrac;
	double		fbar;
	float		mode, median, lquart, uquart, sigma;
	float		*fsample, flquart;
	int		seed = 1;

	fmin = HUGE_VAL; 		/* find max & min */
	fmax = -HUGE_VAL;
	fbar = 0.0;
	ngood = 0;
	badpix = 0;
	for (i = margin; i < N2 - margin; i++) {							
		for (j = margin; j < N1 - margin; j++)	{
			if (fgood(f[i][j])) {
				fmax = (f[i][j] > fmax ? f[i][j] : fmax);
				fmin = (f[i][j] < fmin ? f[i][j] : fmin);
				ngood++;
				fbar += f[i][j];
			} else {
				badpix++;
			}
		}
	}
	srec->badpix = badpix;
	srec->goodpix = ngood;
	if (ngood) {
		srec->fmin = fmin; 
		srec->fmax = fmax;
		fbar /= (float) ngood;
		srec->fmean = (float) fbar;
	} else {
		srec->fmin = 0.0; 
		srec->fmax = 0.0;
		srec->fmean = 0.0;
		srec->flowerquartile = 0.0;
		srec->fupperquartile = 0.0;
		srec->fmedian = 0.0;
		srec->fmode = 0.0;
		srec->sigma = 0.0;
		srec->samplesize = 0;
		return;
	}
	samplefrac = DEFSAMPLEFRAC;		/* take a random sample */
	samplesize = (int) ceil(samplefrac * N1 * N2);
	samplesize = (samplesize > MAX_SAMPLE_SIZE ? MAX_SAMPLE_SIZE : samplesize);
	if (samplesize > ngood) {
		step = 1;
	} else {
		step = (int) floor(sqrt(ngood / (double) samplesize));
	}
	samplesize = (int) (2 + N1 / step) * (2 + N2 / step);

	fsample = (float *) calloc(samplesize, sizeof(float));
	actualsize = 0;
	for (i = margin; i < N2 - margin; i += step) {							
		for (j = margin; j < N1 - margin; j += step)	{
			if (fgood(f[i][j])) {
				fsample[actualsize++] = f[i][j];
			}
			if (actualsize >= samplesize)
				break;
		}
		if (actualsize >= samplesize)
			break;
	}

	srec->samplesize = samplesize = actualsize;

	if (liststats(fsample, samplesize, &median, &lquart, &uquart, &sigma))
		fprintf(stderr, "fdo_stats: warning: problem with stats from sample\n");
	findmode(f, N1, N2, lquart, uquart, &mode, &flquart);
	srec->flowerquartile = lquart;
	srec->fupperquartile = uquart;
	srec->fmedian = median;
	srec->fmode = mode;
	/* empirical fudge */
	srec->sigma = 1.49 * (mode - flquart);
	free(fsample);
}


int			liststats(	float 	*fsample, 
						int		samplesize, 
						float	*medianptr, 
						float	*lquartptr, 
						float	*uquartptr, 
						float	*sigmaptr)
{
	int 	floatcmp();
	int		iplus, iminus, uppi, lowi, imod;
	float	mode, lquart, uquart, median, sigma;
	int	returnvalue = 0;
	
	lowi = floor(0.5 + 0.25 * samplesize);
	uppi = floor(0.5 + 0.75 * samplesize);
	qsort(fsample, samplesize, sizeof(float), floatcmp);
	*lquartptr = lquart = fsample[lowi];
	if (samplesize % 2) {
		*medianptr = median = fsample[samplesize / 2];
	} else {
		*medianptr = median = 0.5 * (fsample[samplesize / 2] + fsample[samplesize / 2 - 1]);
	}
	*uquartptr = uquart = fsample[uppi];

	/* crude estimate of sigma from quartiles */
	*sigmaptr = (uquart - lquart) / (2 * 0.67);
	return (returnvalue);
}




int		floatcmp(float *f1, float *f2)
{
	return (*f1 > *f2 ? 1 : (*f1 == *f2 ? 0 : -1));	/* ascending order */
}

#define	BINS 		1000
#define SIG_FACTOR 	1.0
#define BIG_FACTOR	4

/*
 * - estimate the mode given sorted array of values f with median, quartiles.
 * - works by making a histogram of counts over range BIG_FACTOR * central
 *	50 percential range. 
 *	Smooths with gaussian width ca SIG_FACTOR * sigma (where sigma is
 *	estimated from the quartiles.  
 *
 * we also return 'flquart' which is the vaue such that equally many
 * pixels lie in f < flquart and flquart < f < mode, which is used to give a better sigma..
 *
 */ 
int		findmode(float **f, int N1, int N2, float lquart, float uquart, float *mode, float *flquart)
{
	float 	fmin, fmax;
	long	count[BINS]; 
	int		index, i, j, b, db, bmax, dbmax, submodecount = 0, tempcount = 0;
	float	smoothcount[BINS], smoothcountmax, df, fs, crudesigma;
	float	*expon;
	int	returnvalue;

	if (uquart <= lquart) {
		*mode = 0.5 * (uquart + lquart);
		return (1);
	}
	crudesigma = (uquart - lquart) / (2 * 0.67);
	fs = SIG_FACTOR * crudesigma;
	fmin = 0.5 * (uquart + lquart) - BIG_FACTOR * crudesigma;
	fmax = 0.5 * (uquart + lquart) + BIG_FACTOR * crudesigma;
	df = (fmax - fmin) / BINS;
	dbmax = 3 * fs / df;	/* smoothing window extends to 3-sigma */
	
	expon = (float *) calloc(2 * dbmax + 1, sizeof(float)) + dbmax;
	for (db = - dbmax ; db <= dbmax; db++) {
		expon[db] = exp (-0.5 * df * df * db * db / (fs * fs));
	}	
	
	for (b = 0; b < BINS; b++) 		/* initialise arrays */
		smoothcount[b] = count[b] = 0;
	
	/* make counts */
	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {			
			index = floor(0.5 + (f[i][j] - fmin) / df);
			if (index >= 0 && index < BINS)
				count[index]++;
		}
	}

	/* smooth counts */
	for (b = 0; b < BINS; b++) {
		for (db = - dbmax ; db <= dbmax; db++) {
			if (b + db >= 0 && b + db < BINS)
				smoothcount[b] += count[b + db] * expon[db];
		}	
	}

	smoothcountmax = bmax = 0;				/* find peak of smoothed counts */
	for (b = 0; b < BINS; b++) {
		if (smoothcount[b] > smoothcountmax) {
			smoothcountmax = smoothcount[b];
			bmax = b;
		}
	}
	if (bmax > 0 && bmax < BINS - 1) {		/* seems to have worked */
		*mode = fmin + df * bmax;
		/* compute flquart */
		for (b = 0; b < bmax; b++) {
			submodecount += count[b];
		}
		submodecount += count[bmax] / 2;
		for (b = 0; b < bmax; b++) {
			tempcount += count[b];
			if (tempcount > 0.5 * submodecount)
				break;
		}
		*flquart = fmin + df * b;
		returnvalue = 0;
	} else {
		*mode = 0.5 * (uquart + lquart);
		*flquart = lquart;
		returnvalue = 1;
	}

	

	free(expon - dbmax);
	return (returnvalue);
}


#undef BINS
#undef SIG_FACTOR
#undef BIG_FACTOR
#undef DEFSAMPLESIZE
#undef HUGE_VAL



