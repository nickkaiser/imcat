/*
 * hfindpeaks.c
 *
 * a hierarchical peak finding algorithm
 *
 * creates a sequence of progressively more coarsely smoothed images, finds
 * the peaks at each level of smoothing and links these peaks together into
 * lines according to criteria in neighbours.
 *
 * for each level
 *	* smooth
 *	* find peaks and load them into a list headed by peakhead
 *	  try to install the peaks one at a time onto one of lines in linked list
 *	  headed by linehead.  If we can't we create a new line and link it
 *	  onto the end of the linked list
 *	* find dead lines (ones which didn't get peak added to at last step)
 *	  and dispose of them; i.e. write out most significant point
 */

#define usage "\n\n\n\
NAME\n\
	hfindpeaks --- hierarchical object finder\n\
\n\
SYNOPSIS\n\
	hfindpeaks fitsfile [option...] > catalogue \n\
 		-r rf1 rf2	# range of filter radii (0.5 20.0)\n\
		-d dlnrf	# step in log(rf) (0.2)\n\
		-l flink	# linking parameter(1.0)\n\
		-n nu		# significance threshold (4.0)\n\
		-s sigma mode	# sky statistics\n\
		-a noiseacf	# supply noise autocorrelation function\n\
		-N N1 N2	# set working image size\n\
\n\
DESCRIPTION\n\
	\"hfindpeaks\" --- a heirarchical object finder\n\
	fft gaussian filters an image with sequence of progressively\n\
	larger smoothing radius mexican hat filters and computes\n\
	the significance nu(x; rf) (defined to be the smoothed field at x divided\n\
	by the rms noise fluctuation for smoothing radius rf).\n\
	It then finds peaks of the nu field links these together to\n\
	construct peak trajectories x_pk(rf). We define an object to be the\n\
	point of highest significance along such a trajectory.\n\
	rf1, rf2 are min and max filter radii and we filter\n\
	with logarithmic steps in rf defined by dlnrf\n\
	Peaks at adjacent smoothing levels are connected if their\n\
	separation is less than flink * rf.\n\
	Use -s option to supply sky mode, sigma rather than have\n\
	them calculated from the image.\n\
	The catalogue is in 'lc' format.\n\
	The position 'x' is measured relative to the bottom left\n\
	corner of the bottom left pixel (so e.g. a single 'hot' pixel\n\
	at (ix,iy) = (23,67), would generate an object with x = (23.5, 67.5)\n\
	The rms noise at smoothing scale rf is computed assuming that the\n\
	noise is incoherent.  If the noise is correlated (from resampling\n\
	say) you can supply a noise psf with the -p option.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include "../../fftlib/myfft.h"
#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../utils/args.h"
#include "../../imlib/filters.h"
#include "../../utils/arrays.h"
#include "../../utils/stats_stuff.h"
#include "hfindpeaks.h"

static float	sigma1, sigma2;

static	int	nlines = 0, npeaks = 0;

#define TINY	1.e-20

main(int argc, char *argv[])
{
	float		**acf, **acfs, **fin , **fs, *rf, rf1, rf2, dlnrf, flink, dmax, nulimit, *sigma;
	int		level, nlevels;
	int		icell, jcell, ncells, installed;
	peak		*peakhead, *thepeak, *nextpeak;
	line		*linehead = NULL, *theline, *nextline;
	line		***cell;
	fitsheader	*fits;
	fitscomment	*com;
	char		*fitsfilename, sysstring[1024], tempstring[128], *flag, *noiseacffilename,
			argsstring[128];
	int		kx, ky, i, j, N1, N2, Nmax, M1, M2, NN1, NN2, xoff, yoff;
	FILE		*fitsfile, *noiseacffile;
	float		magicsubstitute;
	fstatsrec	srec;
	float		rfratio, rfactor, ffactor;
	fft_type	fk, fkcopy, acfk;
	int		calcstats, usenoiseacf;
	
	/* defaults */
	rf1 = 0.5;
	rf2 = 20.0;
	dlnrf = 0.2;
	flink = 1.0;
	nulimit = 4.0;
	rfratio = 2.0;		/* mexican hat parameter */
	calcstats = 1;
	usenoiseacf = 0;
	N1 = N2 = 0;

	/* parse args */
	/* intiialise */
	argsinit(argc, argv, usage);
	/* first argument must be a file */
	if (nextargtype() == FLAG_ARG || nextargtype() == NO_ARG) {
		error_exit(usage);
	}
	fitsfilename = getargs();
	if (strcmp(fitsfilename, "-")) {
		fitsfile = fopen(fitsfilename, "r");
		if (!fitsfile) {
			error_exit("hfindpeaks: can't open fits file\n");
		}
	} else {
		fitsfile = stdin;
	}

	/* parse the options */
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'r':
				rf1 = getargf();
				rf2 = getargf();
				break;
			case 'd':
				dlnrf = getargf();
				break;
			case 'l':
				flink = getargf();
				break;
			case 'n':
				nulimit = getargf();
				break;
			case 's':
				srec.sigma = getargf();
				srec.fmode = getargf();
				calcstats = 0;
				break;
			case 'a':
				noiseacffilename = getargs();
				usenoiseacf = 1;
				break;
			case 'N':	
				N1 = getargi();
				N2 = getargi();
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}

	fits = readfitsheader(fitsfile);
	NN1 = fits->n[0];
	NN2 = fits->n[1];
	if (N1) {
		if (NN1 > N1 || NN2 > N2) {
			error_exit("hfindpeaks: image size exceed that specified with -N option\n");
		} 
	} else {
		N1 = NN1;
		N2 = NN2;
	}
	allocFloatArray(&fin, N1, N2);
	for (i = 0; i < N2; i++) {
		for (j = 0; j < N1; j++) {
			fin[i][j] = FLOAT_MAGIC;
		}
	}
	for (i = 0; i < NN2; i++) {
		readfitsline(fin[i], fits);
	}
	close (fitsfile);
	Nmax = (N1 > N2 ? N1 : N2);

	/* allocate the smoothed image */
	allocFloatArray(&fs, N1, N2);

	/* get the mode and sigma if necessary */
	if (calcstats) {
		fdo_stats(fin, NN1, NN2, 0, &srec);
	}
	magicsubstitute = srec.fmode;

	/* create the cathead with lc */
	argsToString(argc, argv, argsstring);
	sprintf(sysstring, "lc -C -b -x -a 'history: %s' -H 'fits_name = {%s}' ", argsstring, fitsfilename);
	sprintf(tempstring, "-H 'fits_size = %d %d 2 vector' -H 'has_sky = 0' ", N1, N2);
	strcat(sysstring, tempstring);
	sprintf(tempstring, "-N '1 2 x' -N '1 1 lg' -N '1 1 rg' -N '1 2 eg' -N '1 1 fs' -N '1 1 nu' < /dev/null");
	strcat(sysstring, tempstring);
	system(sysstring);

	/* replace MAGIC values with mode */
	substitute(fin, N1, N2, magicsubstitute);

	/* compute the rf array */
	nlevels = (int) ceil(log(rf2 / rf1) / dlnrf);
	rf = (float *) calloc(nlevels, sizeof(float));
	for (level = 0; level < nlevels; level++) {
		rf[level] = rf1 * exp(level * dlnrf);
	}

	/* compute the sigma array */

	/* allocate sigma */
	sigma = (float *) calloc(nlevels, sizeof(float));
	/* we do this using fft and images of size 2^N and > 4 * rf2 and >= 32 */
	M1 = 32;
	while (M1 < (int) ceil(4 * rf2)) {
		M1 *= 2;
	}
	M2 = M1;
	allocFloatArray(&acf, M1, M2);
	if (usenoiseacf) {
		/* look at the acf header */
		noiseacffile = fopen(noiseacffilename, "r");
		if (!noiseacffile) {
			error_exit("hfindpeaks: failed to open noise acf image\n");
		}
		fits = readfitsheader(noiseacffile);
		close (noiseacffile);
		/* reread it using makesubimage */
		xoff = (fits->n[0] - M1) / 2;
		yoff = (fits->n[1] - M2) / 2;
		sprintf(sysstring, "makesubimage %d %d %d %d -o < %s", xoff, yoff, M1, M2, noiseacffilename);
		noiseacffile = popen(sysstring, "r");
		if (!noiseacffile) {
			error_exit("hfindpeaks: failed to open noise acf image through makesubimage\n");
		}
		fits = readfitsheader(noiseacffile);
		readfitsplane((void **) acf, fits);
		pclose(noiseacffile);
	} else {
		acf[M2/2][M1/2] = 1;
	}
	allocFloatArray(&acfs, M1, M2);
	alloc_fft(&acfk, M1, M2);
	for (level = 0; level < nlevels; level++) {
		forward_fft(acf, M1, M2, acfk);
		sigma1 = 0.5 * rf[level] * rf[level];
		sigma2 = rfratio * rfratio * sigma1;
		filter(acfk, M1, M2, thefilterfunction);
		filter(acfk, M1, M2, thefilterfunction);
		inverse_fft(acfk, M1, M2, acfs);
		sigma[level] = srec.sigma * sqrt(acfs[M2/2][M1/2] / acf[M2/2][M1/2]);
	}
	

	/* allocate space for fk and a copy of it */
	alloc_fft(&fk, N1, N2);
	alloc_fft(&fkcopy, N1, N2);

	/* make the unsmoothed image fft */
	forward_fft(fin, N1, N2, fk);

	for (level = 0; level < nlevels; level++) {
		dmax = rf[level] * flink;
		ncells = ceil(Nmax / dmax);
		cell = makecell(linehead, dmax, ncells);
		sigma1 = 0.5 * rf[level] * rf[level];
		sigma2 = rfratio * rfratio * sigma1;
		copy_fft(fk, N1, N2, fkcopy);
		filter(fkcopy, N1, N2, thefilterfunction);
		inverse_fft(fkcopy, N1, N2, fs);
		peakhead = getpeaks(N1, N2, fs, level, rf[level], sigma[level]);
		thepeak = peakhead;
		while (thepeak) {			/* loop over new peaks */
			nextpeak = thepeak->next;	/* 'cos thepeak->next gets trashed */
			getcoords(thepeak, &icell, &jcell, dmax);
			installed = 0;
			for (i = icell - 1; i <= icell + 1; i++) {
				for (j = jcell - 1; j <= jcell + 1; j++) {
					if (i < 0 || i >= ncells || j < 0 || j >= ncells)
						continue;
					if (install(thepeak, cell[i][j], dmax)) {
						installed = 1;
						break;
					}
				}
				if (installed)
					break;
			}
			if (!installed)
				addnewline(&linehead, thepeak);
			thepeak = nextpeak;
		}
		removedead(&linehead, level, nulimit);
		freecell(cell, ncells);
		level++;
	}
	
	/* dispose of any lines still alive at end */
	theline = linehead;
	while (theline) {
		nextline = theline->next;
		disposeof(theline, nulimit);
		theline = nextline;
	}
	exit(0);
}



peak	*getpeaks(int N1, int N2, float **fs, int level, float rf, float sigma)
{
	int	i, j;
	peak	*thepeak, *newpeak, *peakhead = NULL;
	float	**frow, **ftop, **fbot, *fp, f, Q11, Q22, Q12, Fp, Fpp;
	
	for (i = 1; i < N2 - 1; i++) {
		frow = fs + i;
		ftop = frow + 1;
		fbot = frow - 1;
		for (j = 1; j < N1 - 1; j++) {
			fp = *frow + j;
			f = *fp;
			/* middle row */
			if (*(fp-1) >= f)
				continue;
			if (*(fp+1) >= f)
				continue;	
			/* top row */
			fp = *ftop + j;
			if (*fp >= f)
				continue;
			if (*(fp-1) >= f)
				continue;
			if (*(fp+1) >= f)
				continue;	
			/* bottom row */
			fp = *fbot + j;
			if (*fp >= f)
				continue;
			if (*(fp-1) >= f)
				continue;
			if (*(fp+1) >= f)
				continue;	
			/* if it survived it is a peak */
			newpeak = (peak *) calloc(1, sizeof(peak));
			npeaks++;
			newpeak->i = i;
			newpeak->j = j;
			newpeak->level = level;
			newpeak->rf = rf;
			newpeak->fs = fs[i][j];
			newpeak->nu = fs[i][j] / sigma;
			Q11 = 2 * fs[i][j] - fs[i+1][j] - fs[i-1][j];
			Q22 = 2 * fs[i][j] - fs[i][j+1] - fs[i][j-1];
			Q12 = 0.25 * (fs[i+1][j+1] - fs[i+1][j-1] - fs[i-1][j+1] + fs[i-1][j-1]);
			newpeak->e1 = - (Q11 - Q22) / (Q11 + Q22);
			newpeak->e2 = 2 * Q12 / (Q11 + Q22);
			/* now we add the simple cludge to get more precise position */
			newpeak->x[0] = newpeak->j + 0.5;
			newpeak->x[1] = newpeak->i + 0.5;
			Fp  = 0.5 * (fs[i][j+1] - fs[i][j-1]);
			Fpp = fs[i][j+1] - 2 * fs[i][j] + fs[i][j-1];
			if (fabs(Fpp) > TINY)
				newpeak->x[0] -= Fp / Fpp;
			Fp  = 0.5 * (fs[i+1][j] - fs[i-1][j]);
			Fpp = fs[i+1][j] - 2 * fs[i][j] + fs[i-1][j];
			if (fabs(Fpp) > TINY)
				newpeak->x[1] -= Fp / Fpp;
			if (!peakhead)
				peakhead = thepeak = newpeak;
			else
				thepeak = thepeak->next = newpeak;
		}
	}
	return (peakhead);
}



int	install(peak *thepeak, line *linehead, float dmax)
{
	line	*theline;
	int	installed;
	
	installed = 0;
	theline = linehead;
	while (theline) {
		if (neighbours(theline->tail, thepeak, dmax)) {	/* add thepeak to the list */
			thepeak->next = NULL;
			theline->tail = theline->tail->next = thepeak;
			installed = 1;
			break;
		}
		theline = theline->cellmate;
	}
	return (installed);
}




void	addnewline(line **lineheadptr, peak *thepeak)
{
	line	*newline;

	newline = (line *) calloc(1, sizeof(line));
	newline->head = newline->tail = thepeak;
	thepeak->next = NULL;
	if (*lineheadptr)
		(*lineheadptr)->prev = newline;
	newline->next = *lineheadptr;
	*lineheadptr = newline;
	nlines++;
}



int	neighbours(peak *peak1, peak *peak2, float dmax)
{
	double	di, dj, dd;
	
	di = peak1->i - peak2->i;
	dj = peak1->j - peak2->j;
	dd = di * di + dj * dj;
	return (dd <= dmax * dmax);
}



void	removedead(line **lineheadptr, int level, float nulimit)
{
	line	*theline, *nextline;
	
	theline = *lineheadptr;
	while (theline) {
		nextline = theline->next;
		if (theline->tail->level != level) {	/* a dead one */
			if (theline == *lineheadptr) {
				*lineheadptr = nextline;
			} else {
				theline->prev->next = nextline;
			}
			if (nextline)
				nextline->prev = theline->prev;
			disposeof(theline, nulimit);
		}
		theline = nextline;
	}
}



void	disposeof(line *theline, float nulimit)
{
	peak	*thepeak, *nextpeak, *bestpeak;
	float	maxsig = 0.0;
	
	/* find the highest peak */
	thepeak = theline->head;
	while (thepeak) {
		nextpeak = thepeak->next;
		if (thepeak->nu > maxsig) {
			maxsig = thepeak->nu;
			bestpeak = thepeak;
		}
		thepeak = nextpeak;
	}

	/* write the best peak if it qualifies */
	if (maxsig > nulimit) {
		thepeak = theline->head;
		while (thepeak) {
			if (thepeak->next == bestpeak && bestpeak->next != NULL) {
				output(bestpeak, thepeak, bestpeak->next);
				break;
			}
			thepeak = thepeak->next;
		}
	}

	/* free up memory */
	thepeak = theline->head;
	while (thepeak) {
		nextpeak = thepeak->next;
		free(thepeak);
		npeaks--;
		thepeak = nextpeak;
	}
	free(theline);
	nlines--;
}




line	***makecell(line *linehead, float dmax, int ncells)
{
	line 	***cell, *theline;
	int	icell, jcell;
		
	cell = (line ***) calloc(ncells, sizeof(line **));
	for (icell = 0; icell < ncells; icell++)
		cell[icell] = (line **) calloc(ncells, sizeof(line *));
	theline = linehead;
	while (theline) {
		getcoords(theline->tail, &icell, &jcell, dmax);
		if (icell < 0 || icell >= ncells || jcell < 0 || jcell >= ncells) {
			fprintf(stderr, "whoops\n");
			exit(-1);
		}
		theline->cellmate = cell[icell][jcell];
		cell[icell][jcell] = theline;
		theline = theline->next;
	}
	return (cell);
}



void	freecell(line ***cell, int ncells)
{
	int icell;
	for (icell = 0; icell < ncells; icell++) {
		free(cell[icell]);
	}
	free(cell);
}


void	getcoords(peak *thepeak, int *icell, int *jcell, float dmax)
{
	*icell = floor(thepeak->i / dmax);
	*jcell = floor(thepeak->j / dmax);
}


#define ZMAX	10

/*
 * mexican hat filter function 
 */
float	thefilterfunction(float ki, float kj)
{
	float	kk, z1, z2;
		
	kk = ki * ki + kj * kj;
	z1 = kk * sigma1;
	if (z1 > ZMAX)
		return (0.0);
	z2 = kk * sigma2;
	if (z2 > ZMAX) {
		return (exp(-z1));
	} else {
		return (exp(-z1) - exp(-z2));
	}
}

#undef ZMAX


float	interp(float d, float x1, float x2, float x3)
{
	if (d > 0)
		return((1 - d) * x2 + d * x3);
	else
		return((1 + d) * x2 - d * x1);
}

/* empirical factors to convert rg => rh, fs => mag */
#define RFACTOR         0.66
#define FFACTOR         15.41
#define EFACTOR         1.50
#define	BUFFSIZE	8

void	output(peak *bestpk, peak *prevpk, peak *nextpk)
{
	float		d, nu1, nu2, nu3;
	static double	opbuff[BUFFSIZE];

	nu1 = prevpk->nu;
	nu2 = bestpk->nu;
	nu3 = nextpk->nu;
	d = 0.5 * (nu3 - nu1) / (2 * nu2 - nu1 - nu3);
	d = (d > 1 ? 1 : d);
	d = (d < -1 ? -1 : d);
	opbuff[0] = (double) bestpk->x[0];
	opbuff[1] = (double) bestpk->x[1];
	opbuff[2] = (double) (FFACTOR * bestpk->fs * bestpk->rf * bestpk->rf);
	opbuff[3] = (double) (RFACTOR * interp(d, prevpk->rf, bestpk->rf, nextpk->rf));
	opbuff[4] = (double) (-1.0 * EFACTOR * bestpk->e1);
	opbuff[5] = (double) (EFACTOR * bestpk->e2);
	opbuff[6] = (double) interp(d, prevpk->fs, bestpk->fs, nextpk->fs);
	opbuff[7] = (double) interp(d, prevpk->nu, bestpk->nu, nextpk->nu);
	fwrite(opbuff, sizeof(double), BUFFSIZE, stdout);
}

#undef	BUFFSIZE
#undef	RFACTOR
#undef	FFACTOR
#undef	EFACTOR
