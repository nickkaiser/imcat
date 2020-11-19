#define	usage "\n\n\n\
SYNOPSIS\n\
		analyse	[option...] < input.cat > output.cat\n\
			-n nu		# threshold nu (default -10.0)\n\
			-a a1 a2	# annulus for sky: a1 < r <= a2, (default 16 32)\n\
			-f		# fit objects to Gaussian ellipsoids\n\
			-R r0		# gaussian window scale (3.0)\n\
			-x mul		# gaussian window scale  = mul * obj->rg (1.0)\n\
			-e n rc		# use W(r) = (r^2 + rc^2)^n/2 window (-2, 2)\n\
			-s		# force recalculation of local sky params\n\
			-z		# switch off local sky determination\n\
			-r alpha	# aperture = alpha * r_petrosian (3.0)\n\
			-m x a		# zap neighbouring images\n\
			-S		# force recalculation of image sigma, mode\n\
			-Q		# sky annulus matched to aperture\n\
			-F deltam	# filter output to reject bad rh, mag\n\
\n\
DESCRIPTION\n\
		\"analyse\" analyses images around a catalogue of objects\n\
		created by (h)findpeaks. It determines a constant plus gradient model\n\
		for the local sky parameters using NE, NW, SW, SE quadrants of\n\
		an annulus (unless you tell it not to).\n\
		Only uses objects with nu (determined by findpeaks) above threshold.\n\
		-f to do gaussian ellipsoid fit.\n\
		By default, moments etc are determined with gaussian window scale 3.0 pixels.\n\
		-e to specify a softened power law window.\n\
		-R option to change scale length.\n\
		-x option to override this and set window scale for moments\n\
		to be mul * obj.rg (rg as determined by hcat2cat)\n\
		Use -ve alpha to use aperture = (-alpha) * obj->rg for photometry.\n\
		Luminosities incorporate 1 / normfactor from fits header\n\
		-m x a to zap disks radius r * r_x, where x = \"g\" or \"n\"\n\
		(for r_numax), around neighbours\n\
		-Q with -ve alpha to make reference annlus run from r_ap to 2 * r_ap\n\
		With -F option, if total mag differs from hfindpeaks mag estimate by > deltam\n\
		we ignore rh, total mag in favor of hfindpeaks values.\n\
\n\n\n"		

/* modified do_object_stats, formerly in object_stuff.c 
July 3, 1992:
Current version works as follows:
July 3
	* calculate the medians, non-MAGIC count and centroids for 4 surronding 
	  sectors (NSEW) in an annulus a1 < r <= a2.   
	* for each peak, accumulate 1 and f in unit width rings 
		out to max radius GC_MAX (GC -> "growth curve").
	* find r_numax = integer radius of max sigma for interior light
	* set rmax = alpha * r_numax (or GC_MAX - 1, whichever is smaller)
	* calculate total luminosity, halflight radius, various moments and e1, e2
Mon Jul  6
	* added mode estimator for surrounding sectors
Wed Aug  5 23:38:03 EDT 1992
	* power law window function for ellipticity
Thu Aug  6 10:05:40 EDT 1992
	* new cat version - comments, hassky etc.
Fri Dec 10 09:37:13 EST 1993
	added estimators for Psm11, Psm22 and Psh = (Psh11 + Psh22) / 2
Sat Dec 25 11:17:22 EST 1993
	changed to gaussian window for moments
Sun Sep  4 16:19:14 EDT 1994
	zapping neighbours added
		this works by creating new images fzap = f and nzap = 0:
			for each object
			do
				for each pixel in disk
				do
					set fzap = MAGIC
					increment nzap by 1
				done
			done
		we can then restoreobj() unmasked parts of a target object to fzap
		by copying disk pixels from f, but only if nzap = 1 and then
		eraseobj() to reset those pixels to MAGIC
		we then analyse fzap.
	converted to a pipe
Sat Nov 19 13:24:05 EST 1994
	converted to internal float image format
*/
	

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "fits.h"
#include "error.h"
#include "arrays.h"
#include "stats_stuff.h"
#include "object_stuff.h"
#include "analyse.h"
#include "magic.h"
#include "catio.h"
#include "gaussfit/gaussfit.h"
#include "status.h"


#define RMAX_MIN	6

float		ff, ffi, ffj, ffmode;			/* parameters for fsky() model */
int			fit = 0;
int		powerlawwindow;
int		nosky = 0;

#define G_RADIUS 0
#define N_RADIUS 1

#define	ZAP	0
#define RESTORE	1

main(int argc, char *argv[])	
{
	float 		**fzap, **f;
	short 		**nzap;
	int		arg = 1, N1, N2, a1, a2, s, i, j; 
	fstatsrec	srec;
	object 		obj;
	cathead		thecat;
	FILE		*fitsf, *tempf;
	int		com, comc;
	char		*comv[MAX_COMMENTS];
	float		nu;
	skyquad		sky;
	char		systemstring[1024];
	float		ne, alpha, rc, r0, azap;
	int		alreadyhassky, redosky, objnumber = 0, scaledwindow, Ffilter;
	int		zapneighbours, zapradiustype, redostats, margin, matchedannulus;
	float		normfactor, mul, deltam;
	
	/* defaults */
	a1 = 16;
	a2 = 32;
	nu = -10.0;
	ne = -2.0;
	redosky = 0;
	alpha = 3.0;
	rc = 2;
	powerlawwindow = 0;
	r0 = 3.0;
	mul = 1.0;
	scaledwindow = 0;
	zapneighbours = 0;
	redostats = 0;
	margin = 3;
	matchedannulus = 0;
	Ffilter = 0;
	
	while (arg < argc) {
		if (*argv[arg] == '-') {
			switch (*(argv[arg++]+1)) {
				case 'e':
					if (EOF == sscanf(argv[arg++], "%f", &ne))
						error_exit(usage);
					if (EOF == sscanf(argv[arg++], "%f", &rc))
						error_exit(usage);
					powerlawwindow = 1;
					if (scaledwindow)
						error_exit("analyse: cannot use -x and -e options together\n");
					break;
				case 'n':
					if (EOF == sscanf(argv[arg++], "%f", &nu))
						error_exit(usage);
					break;
				case 'r':
					if (EOF == sscanf(argv[arg++], "%f", &alpha))
						error_exit(usage);
					break;
				case 'R':
					if (EOF == sscanf(argv[arg++], "%f", &r0))
						error_exit(usage);
					break;
				case 'a':
					if (EOF == sscanf(argv[arg++], "%d", &a1))
						error_exit(usage);
					if (EOF == sscanf(argv[arg++], "%d", &a2))
						error_exit(usage);
					redosky = 1;
					break;
				case 'f':
					fit = 1;
					break;
				case 'z':
					nosky = 1;
					break;
				case 's':
					redosky = 1;
					break;
				case 'S':
					redostats = 1;
					break;
				case 'Q':
					matchedannulus = 1;
					break;
				case 'F':
					Ffilter = 1;
					if (EOF == sscanf(argv[arg++], "%f", &deltam))
						error_exit(usage);
					break;
				case 'm':
					zapneighbours = 1;
					switch (argv[arg++][0]) {
						case 'g':
							zapradiustype = G_RADIUS;
							break;
						case 'n':
							zapradiustype = N_RADIUS;
							break;
						default:
							fprintf(stderr, "analyse: bad zap radius type\n");
							error_exit(usage);
							break;
					}
					if (EOF == sscanf(argv[arg++], "%f", &azap))
						error_exit(usage);
					break;
				case 'x':
					if (EOF == sscanf(argv[arg++], "%f", &mul))
						error_exit(usage);
					scaledwindow = 1;
					if (powerlawwindow)
						error_exit("analyse: cannot use -x and -e options together\n");
					break;
				default:
					error_exit(usage);
					break;
			}
		} else {
			error_exit(usage);
		}
	}
	if (alpha >= 0.0 && matchedannulus)
		error_exit("analyse: use negative alpha with -Q option\n");

	read_cat_head(&thecat);
	alreadyhassky = thecat.hassky;
	thecat.hassky = 1;
	fitsf = fopen(thecat.image, "r");
	if (!fitsf)
		error_exit("analyse: unable to open fits file for input");
	set_fits_ipf(fitsf);
	fread_fits(&f, &N1, &N2, &comc, comv);
	/* get normfactor */
	normfactor = 1.0;
	for (com = 0; com < comc; com++) {
		if (0 == strncmp("HISTORY = normalise fac=", comv[com], 23)) {
			sscanf(comv[com], "HISTORY = normalise fac= %f", &normfactor);
		}
	}
	if (redostats)
		fdo_stats(f, N1, N2, margin, &(thecat.srec));
	srec = thecat.srec;
	ffmode = srec.fmode;
	fclose(fitsf);

	/* to zap neighbours we need to do a first pass to make zapped images */
	if (zapneighbours) {
		allocFloatArray(&fzap, N1, N2);
		for (i = 0; i < N2; i++)
			for (j = 0; j < N2; j++)
				fzap[i][j] = f[i][j];
		allocShortArray(&nzap, N1, N2);
		tempf = fopen("analyse.tmp", "w");
		if (!tempf)
			error_exit("analyse: fopen() failed\n");
		set_cat_opf(tempf);
		write_cat_head(&thecat);
		while (read_object(&obj)) {
			zap(ZAP, &obj, zapradiustype, azap, f, fzap, nzap, N1, N2);
			write_object(&obj);
		}
		fclose(tempf);
		tempf = fopen("analyse.tmp", "r");
		if (!tempf)
			error_exit("analyse: fopen() failed\n");
		set_cat_ipf(tempf);
		read_cat_head(&thecat);
		set_cat_opf(stdout);
	} else {
		fzap = f;
	}

	add_cat_comment(argc, argv, &thecat);
	write_cat_head(&thecat);
	while (read_object(&obj)) {
		if (obj.nu > nu) {
			if (zapneighbours)
				zap(RESTORE, &obj, zapradiustype, azap, f, fzap, nzap, N1, N2);
			if (!alreadyhassky || redosky) {
				if (matchedannulus) {
					a1 = ceil(-alpha * obj.rg);
					a2 = 2 * a1;
				}
				dosky(fzap, N1, N2, obj.i, obj.j, a1, a2, &sky, srec.fmode, srec.sigma);
				obj.sky = sky;
			} else {
				sky = obj.sky;
			}
			setskyparameters(&sky);
			if (powerlawwindow)
				do_object_stats(&obj, fzap, N1, N2, fsky, srec.sigma, ne, rc, alpha);
			else
				if (scaledwindow)
					do_object_stats(&obj, fzap, N1, N2, fsky, srec.sigma, 0.0, mul * obj.rg, alpha);
				else
					do_object_stats(&obj, fzap, N1, N2, fsky, srec.sigma, 0.0, r0, alpha);
			if (zapneighbours)
				zap(ZAP, &obj, zapradiustype, azap, f, fzap, nzap, N1, N2);
			/* now correct luminosity */
			obj.l /= normfactor;
			if (Ffilter) {
				if (obj.l > obj.lg * exp(0.92 * deltam) || obj.l < obj.lg * exp(-0.92 * deltam)) {
					obj.l = obj.lg;
					obj.rh = RFACTOR * obj.rg;
				}
			}
			write_object(&obj);
		}
	}
	fclose(tempf);
	exit(0);
}

void	dosky(	float **f, int N1, int N2, int ip, int jp, int a1, int a2, 
	skyquad *sky, float fmode, float sigma)
/**
 ** calculates occupation and mode for NSEW sectors;
 **/
{
	int 					i, j, m, ii, s;
	static int				lasta2 = 0;
	static float 			*fs[4] = {NULL, NULL, NULL, NULL};
	int						sortarraysize;
	float					rdefault;
	float					mode;
	float					median, lquart, uquart, quadsigma;
						
	sortarraysize = (a2 + 1) * (a2 + 1);
	if (a2 != lasta2) {			/* need to allocate new sort arrays */
		lasta2 = a2;
		for (s = 0; s < 4; s++)	{
			if (fs[s] != NULL)	/* free the old ones if they exist */
				free(fs[s]);
			fs[s] = (float *) calloc(sortarraysize, sizeof(float));
			if (!fs[s])
				error_exit("dosky: memory allocation failed");
		}
	}
	
	rdefault = 8 * (a2 * a2 * a2 - a1 * a1 * a1) / (3 * sqrt(2.0) * PI *
		(a2 * a2 - a1 * a1));
	
	for (s = 0; s < 4; s++) {	
		sky->i[s] = sky->j[s] = sky->n[s] = 0;		/* zero cumulants */
		for (m = 0; m < sortarraysize; m++) 		/* sero sort arrays */
			fs[s][m] = 0;
	}
	
	for (i = ip - a2; i <= ip + a2; i++) {
		if (i < 0 || i >= N2)
			continue;
		for (j = jp - a2; j <= jp + a2; j++) {
			if (j < 0 || j >= N1)
				continue;
			ii = (i - ip) * (i - ip) + (j - jp) * (j - jp);
			if (ii > a2 * a2 || ii <= a1 * a1)
				continue;
			if (f[i][j] != MAGIC) {
				if ((j - jp) > 0 && (i - ip) >= -(j - jp) && (i - ip) < (j - jp))
					s = NORTH;
				if ((j - jp) < 0 && (i - ip) > (j - jp) && (i - ip) <= -(j - jp)) 			
					s = SOUTH;
				if ((i - ip) > 0 && (j - jp) > -(i - ip) && (j - jp) <= (i - ip))
						s = EAST;
				if ((i - ip) < 0 && (j - jp) >= (i - ip) && (j - jp) < -(i - ip))
						s = WEST;
				if (sky->n[s] >= sortarraysize)
						error_exit("dosky: this cannot happen");
				fs[s][sky->n[s]] = f[i][j];
				sky->n[s]++;
				sky->i[s] += i;
				sky->j[s] += j;
			}	
		}
	}
	
	for (s = 0; s < 4; s++) {		/* normalise centroids */
		if (sky->n[s]) {					
			sky->i[s] /= sky->n[s]; 
			sky->j[s] /= sky->n[s];
			liststats(fs[s], sky->n[s], &mode, &median, &lquart, &uquart, &quadsigma);
			sky->f[s] = mode;
		} else {				/* default behaviour for empty sectors */
			sky->f[s] = fmode;
			switch (s) {
				case NORTH:
					sky->i[s] = 0;
					sky->j[s] = rdefault;
					break;
				case SOUTH:
					sky->i[s] = 0;
					sky->j[s] = - rdefault;
					break;
				case EAST:
					sky->i[s] = rdefault;
					sky->j[s] = 0;
					break;
				case WEST:
					sky->i[s] = - rdefault;
					sky->j[s] = 0;
					break;
				default:
					error_exit("this cannot happen\n");
					break;
			}
		}
	}
}

/*
 * ne is the power law index for the ellipticity window function
 * aperture is alpha * r_numax
 */
void	do_object_stats(object *pk, float **f, int N1, int N2, float(*fsky)(int i, int j), float sigma,
			float ne, float rc, float alpha)
{
	static 	float	sumf[GC_MAX], sum1[GC_MAX], summagic[GC_MAX];
	float		q11, q22, q12, d[2], denom;
	int		r, i, j, di, dj, rmax;
	float		dx, dy, f0, a, b, phi;
	static	float 	**ffit = 0;
	int		gfstatus = 0;
	float		W, Wp, Wpp, fc, DD, DD1, DD2;
	float		Xsm11, Xsm22, Xsm12, Xsh11, Xsh22, Xsh12, em[2], eh[2];	/* smear polarizability bits */
	float		temp;
	
	pk->status = 0;

	for (r = 0; r < GC_MAX; r++)					/* zero the sums */
		sumf[r] = sum1[r] = summagic[r] = 0;
		
	for (i = pk->i - GC_MAX; i <= pk->i + GC_MAX; i++){		/* sum over shells */
		if (i < 0 || i >= N2) continue;
		di = i - pk->i;
		for (j = pk->j - GC_MAX; j <= pk->j + GC_MAX; j++){
			if (j < 0 || j >= N1) continue;
			dj = j - pk->j;
			r = floor(0.5 + sqrt((double) (di * di + dj * dj)));
			if (r >= 0 && r < GC_MAX){
				if (f[i][j] == MAGIC) {
					summagic[r]++;
				} else {
					sum1[r] 	+= 1;
					sumf[r] 	+= (f[i][j] - fsky(di, dj));
				}
			}
		}
	}

	pk->rnumax = GC_MAX;	
	for (r = 1; r < GC_MAX; r++) {				/* cumulate sums */		
		sum1[r] += sum1[r - 1];			/* and find max-nu radius */
		sumf[r] += sumf[r - 1];
		summagic[r] += summagic[r - 1];
		if (sumf[r] * sumf[r] * sum1[r-1] < sumf[r-1] * sumf[r-1] * sum1[r]
				&& pk->rnumax == GC_MAX) {
			pk->rnumax = r - 1;			/* peak-nu radius*/
			if (sum1[r-1] > 0)
				pk->numax = sumf[r - 1] / (sqrt((double) sum1[r - 1]) * sigma);
			else
				pk->numax = 0;
		}
	}
	
	if (alpha > 0) {
		rmax = floor(0.5 + (alpha * pk->rnumax > GC_MAX -1 ? GC_MAX - 1 : alpha * pk->rnumax));
		if (rmax < RMAX_MIN)
			rmax = RMAX_MIN;
	} else {
		rmax = floor(0.5 - alpha * pk->rg);
		rmax = (rmax > GC_MAX -1 ? GC_MAX - 1 : rmax);
	}
	pk->rmax = rmax;
	pk->l = sumf[rmax];
/*	pk->sigmal = sqrt((double) sum1[rmax]) * sigma; discontinued */
	pk->nmagic = summagic[rmax];
		
	for (r = 0; r < GC_MAX; r++)		/* do half light radius */
		if (sumf[r] > 0.5 * pk->l)
			break;
	if (r == 0) {
/*		pk->rh = 0.1; we keep value from hfindpeaks */
		pk->status += TINY_RADIUS;
	} else {
		pk->rh = (r - 1);
		if (sumf[r] - sumf[r-1] > 0.0)
			pk->rh += (0.5 * pk->l - sumf[r-1]) / (sumf[r] - sumf[r-1]);
	}
	
	/* now do 1st and second moments and Psm */
	q11 = q22 = q12 = d[0] = d[1] = pk->invr2 = 0.0;
	Xsm11 = Xsm22 = Xsm12 = Xsh11 = Xsh22 = Xsh12 = em[0] = em[1] = eh[0] = eh[1] = 0.0;
	for (i = pk->i - rmax; i <= pk->i + rmax; i++) { 
		for (j = pk->j - rmax; j <= pk->j + rmax; j++) {
			if (i >= 0 && i < N2 && j >= 0 && j < N1) {
				di = i - pk->i;
				dj = j - pk->j;
				r = floor(0.5 + sqrt((double) (di * di + dj * dj)));
				if (r <= rmax) {
					if (f[i][j] != MAGIC) {
						if (powerlawwindow) {
							W = pow((double) (r * r + rc * rc), (double) (0.5 * ne));
							Wp = 0.5 * ne * pow((double) (r * r + rc * rc), (double) (0.5 * ne - 1));
							Wpp = 0.5 * ne * (0.5 * ne - 1) * 
								pow((double) (r * r + rc * rc), (double) (0.5 * ne - 2));;
						} else {
							W = exp(-0.5 * r * r / (rc * rc));
							Wp = -0.5 * W / (rc * rc);
							Wpp = 0.25 * W / (rc * rc * rc * rc);
						}
						fc = (f[i][j] - fsky(di, dj));
						d[0] += fc * di;
						d[1] += fc * dj;
						q11 += fc * W * di * di;
						q22 += fc * W * dj * dj;
						q12 += fc * W * di * dj;
						DD = di * di + dj * dj;
						DD1 = di * di - dj * dj;
						DD2 = 2 * di * dj;
						Xsm11 += (2 * W + 4 * Wp * DD + 2 * Wpp * DD1 * DD1) * fc; 
						Xsm22 += (2 * W + 4 * Wp * DD + 2 * Wpp * DD2 * DD2) * fc; 
						Xsm12 += 2 * Wpp * DD1 * DD2 * fc; 
						Xsh11 += (2 * W * DD + 2 * Wp * DD1 * DD1) * fc;
						Xsh22 += (2 * W * DD + 2 * Wp * DD2 * DD2) * fc;
						Xsh12 += 2 * Wp * DD1 * DD2 * fc;
						em[0] += (6 * Wp + 2 * Wpp * DD) * DD1 * fc;
						em[1] += (6 * Wp + 2 * Wpp * DD) * DD2 * fc;
						eh[0] += 2 * Wp * DD * DD1 * fc;
						eh[1] += 2 * Wp * DD * DD2 * fc;
						pk->invr2 += (f[i][j] - fsky(di, dj)) * W;
					}
				}
			}
		}
	}
	/* calculate ellipticities */
	denom = q11 + q22;			/* can change to sqrt(det) if desired */
	if (denom > 0) {
		pk->invr2 /= denom;
		pk->e[0] = (q11 - q22) / denom;
		pk->e[1] = 2 * q12 / denom;
		em[0] /= denom;
		em[1] /= denom;
		eh[0] /= denom;
		eh[1] /= denom;
		eh[0] += 4 * pk->e[0];
		eh[1] += 4 * pk->e[1];
		pk->Psm11 = Xsm11 / denom - pk->e[0] * em[0];
		pk->Psm22 = Xsm22 / denom - pk->e[1] * em[1];
		pk->Psm12 = Xsm12 / denom - 0.5 * (pk->e[0] * em[1] + pk->e[1] * em[0]);	/* fudge! */
		pk->Psh11 = Xsh11 / denom - pk->e[0] * eh[0];
		pk->Psh22 = Xsh22 / denom - pk->e[1] * eh[1];
		pk->Psh12 = Xsh12 / denom - 0.5 * (pk->e[0] * eh[1] + pk->e[1] * eh[0]);
	} else {
/*		pk->e[0] = pk->e[1] = 0;*/
		pk->Psm11 = pk->Psm22 = pk->Psm12 = pk->Psh11 = pk->Psh22 = pk->Psh12 = 0.0;
		pk->status += NEG_QTRACE;
		pk->invr2 = 0;
	}

	if (fit) {				/* fitting */
		if (!ffit) {		/* first time around allocate ffit */
			ffit = (float **) calloc(2 * GC_MAX + 1, sizeof(float *)) + GC_MAX;
			for (i = -GC_MAX; i <= GC_MAX; i++)
				ffit[i] = (float *) calloc(2 * GC_MAX + 1, sizeof(float)) + GC_MAX;
		}
		for (di = -GC_MAX; di <= GC_MAX; di++) {	/* make a copy */
			for (dj = -GC_MAX; dj <= GC_MAX; dj++) {
				if (di < -rmax || di > rmax || dj < -rmax || dj > rmax) {
					ffit[di][dj] = 0;
				} else {
					i = pk->i + di;
					j = pk->j + dj;
					if (i < 0 || i >= N2 || j < 0 || j >= N1)
						ffit[di][dj] = 0;
					else
						ffit[di][dj] = f[i][j] - fsky(di, dj);
				}
			}
		}
		/* starting guess
		if (pk->rh > 0.1) {
			f0 = pk->l / (2 * PI * pk->rh * pk->rh);
			a = pk->rh;
		} else { */
			f0 = ffit[0][0];
			a = 1;
	/*	}*/
		gfstatus = gaussfit(ffit, rmax, 0, 0, &dx, &dy, &a, &b, &phi, &f0);
		if (f0 < 0)
			gfstatus += NEG_FLUX;
		if (fabs(dx) > PEAK_WAND_LIM || fabs(dx) > PEAK_WAND_LIM)
			gfstatus += PEAK_WANDERED;
		if (a <= 0 || b <= 0)
			gfstatus += NEG_RADIUS;
		pk->status = gfstatus;
		if (!gfstatus) {			/* fitting was good */
			pk->rh = sqrt(a * b);
			pk->l = f0 * 2 * PI * a * b;
			pk->e[0] = (1 - b / a) * cos(2 * phi);
			pk->e[1] = (1 - b / a) * sin(2 * phi);
		}
	}
}




/* next two functions define the sky model - here a simple straight
average plus linear trend */

void	setskyparameters(skyquad *sky)
{
	ff = 	0.25 * (sky->f[NORTH] + sky->f[SOUTH] + sky->f[EAST] + sky->f[WEST]);
	ffi = 	(sky->f[EAST] - sky->f[WEST]) / (sky->i[EAST] - sky->i[WEST]);
	ffj = 	(sky->f[NORTH] - sky->f[SOUTH]) / (sky->j[NORTH] - sky->j[SOUTH]);
}

float	fsky(int di, int dj)
{
	if (nosky)
		return(ffmode);
	else
		return (ff + di * ffi + dj * ffj);
}




void	zap(int zapmode, object *obj, int radiustype, float a, float **f, float **fzap, short **nzap, int N1, int N2)
{
	/* set fzap to MAGIC in disk and increment nzap */
	int	i, j, d;
	float	rmax, rr, rrmax;
	
	switch (radiustype) {
		case G_RADIUS:
			rmax = a * obj->rg;
			break;
		case N_RADIUS:
			rmax = a * obj->rnumax;
			break;
		default:
			error_exit("zap: bad radiustype\n");
			break;
	}
	d = ceil(rmax);
	rrmax = rmax * rmax;
	for (i = obj->i - d; i <= obj->i + d; i++) {
		if (i < 0 || i >= N2)
			continue;
		for (j = obj->j - d; j <= obj->j + d; j++) {
			if (j < 0 || j >= N1)
				continue;
			rr = (i - obj->i) * (i - obj->i) + (j - obj->j) * (j - obj->j);
			if (rr <= rrmax) {
				switch (zapmode) {
					case ZAP:
						fzap[i][j] = MAGIC;
						nzap[i][j]++;
						break;
					case RESTORE:
						if (nzap[i][j] == 1)
							fzap[i][j] = f[i][j];
						nzap[i][j]--;
						break;
					default:
						error_exit("zap: bad zapmode\n");
						break;
				}
			}
		}
	}
}



