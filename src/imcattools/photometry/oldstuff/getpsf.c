/*
 * getpsf.c
 */


#define	usage "\n\n\n\
SYNOPSIS\n\
		getpsf [options...] a.cat b.cat ......\n\
			-s		# generate getpsf.out\n\
\n\
DESCRIPTION\n\
		\"getpsf\" fits stars in specified range to gaussians\n\
		Invokes sm to popup a window with r-l scatterplot so user\n\
		can define sample of stars - then fills rpsf field in xxx.cat\n\
	        with -s option we analyse -f the stars and psf parameters\n\
		are written to getpsf.out\n\
\n\n\n"		


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "fits.h"
#include "error.h"
#include "stats_stuff.h"
#include "object_stuff.h"
#include "catio.h"
#include "smpopup.h"

#define       MAX(x,y) (((x) > (y)) ? (x) : (y))
#define       MIN(x,y) (((x) < (y)) ? (x) : (y))

#define		MAX_STARS	1000
#define		RLMODE		0
#define		EEMODE		1

#define		ELIMIT 		0.3

float		*smr, *sml, smx1 = 0, smy1 = 0, smx2 = 0, smy2 = 0;
int			smN, smmode, nstars;
float		e1[MAX_STARS], e2[MAX_STARS], l[MAX_STARS], rh[MAX_STARS];
int			i[MAX_STARS], j[MAX_STARS];
void		drawselectionbox();

main(int argc, char *argv[])	
{
	object 		obj;
	cathead		thecat;
	FILE		*catf, *opf, *paramf;
	int		comc;
	char		*comv[MAX_COMMENTS], catfilename[128];
	float		e, a, b, phi;
	int		iobj, nobjects, istar, ngood, ncats, icat;
	float		rmin, rmax, lmin, lmax, e1min, e1max, e2min, e2max, e1av, e2av, rav;
	char		systemstring[1024];
	int		nargs = 1, makescript = 0;
	
	if (argv[1][0] == '-')
		if (argv[1][1] == 's') {
			nargs++;
			makescript = 1;
		} else
			error_exit(usage);

	ncats = argc - nargs;

	if (makescript) {
		paramf = fopen("getpsf.out", "w");
		fprintf(paramf, "#        a          b        phi catfilename\n");
	}

	for (icat = nargs; icat < ncats + nargs; icat++) {
	/* open the cat file */
	strcpy(catfilename, argv[icat]);
	fprintf(stderr, "#\n# processing %s\n", catfilename);
	if (!(catf = fopen(catfilename, "r")))
		error_exit("psfcorrect: failed to open cat file\n");
	
	/* count the objects */
	set_cat_ipf(catf);
	read_cat_head(&thecat);
	nobjects = 0;
	while (read_object(&obj)) {
		nobjects++;
	}
	fclose(catf);
	
	/* allocate and fill the r, l arrays */
	smr = (float *) calloc(nobjects, sizeof(float));
	sml = (float *) calloc(nobjects, sizeof(float));
	catf = fopen(catfilename, "r");
	read_cat_head(&thecat);
	iobj = 0;
	while (read_object(&obj)) {
		smr[iobj] = obj.rh;
		sml[iobj] = (obj.l > 0.0 ? log(obj.l) : 0.0);
		iobj++;
	}
	smN = nobjects;
	fclose(catf);
	
	/* make r-l scatter plot & get user to select stars */
	fprintf(stderr, "#\n# use cursor + key 'b' until you get a nice box around stars\n");
	fprintf(stderr, "# then hit 'x' to exit from cursor input mode\n#\n");
	smmode = RLMODE;
	smpopup(drawfn, cursorfn, "-g 400x400+100+100");
	rmin = MIN(smx1, smx2);
	rmax = MAX(smx1, smx2);
	lmin = MIN(smy1, smy2);
	lmax = MAX(smy1, smy2);
	fprintf(stderr, "# selected range\n#\tl = %10.3e to %10.3e\n#\tr = %10.3e to %10.3e\n",
		lmin, lmax, rmin, rmax);
	sprintf(systemstring, "select -l %.3f %.3f -rh %.3f %.3f < %s > getpsf.tmp",
		lmin, lmax, rmin, rmax, catfilename);
	fprintf(stderr, "# systemstring = %s\n", systemstring);
	if (0 > system(systemstring))
		error_exit("getpsf: system call 1 failed\n");
		
	/* calculate average rpsf */
	if (!(catf = fopen("getpsf.tmp", "r")))
		error_exit("psfcorrect: failed to open getpsf.tmp\n");
	set_cat_ipf(catf);
	read_cat_head(&thecat);
	nstars = 0;
	rav = 0.0;
	while (read_object(&obj)) {
		nstars++;
		rav += obj.rh;
	}
	fclose(catf);
	if (nstars)
		rav /= nstars;
	else
		error_exit("getpsf: you fool, you selected no stars!!\n");
	fprintf(stderr, "# %d stars averaged: <rh> = %.3f\n", nstars, rav);
	
	/* rewrite original cat with rpsf field filled */
	sprintf(systemstring, "/bin/cp %s getpsf.tmp2", catfilename);
	fprintf(stderr, "# systemstring = %s\n", systemstring);
	if (0 > system(systemstring))
		error_exit("getpsf: system call 2 failed\n");
	if (!(catf = fopen("getpsf.tmp2", "r")))
		error_exit("psfcorrect: failed to open getpsf.tmp2\n");
	set_cat_ipf(catf);
	if (!(opf = fopen(catfilename, "w")))
		error_exit("psfcorrect: failed to reopen original cat to enter rpsf\n");
	set_cat_opf(opf);
	read_cat_head(&thecat);
	thecat.rpsf = rav;
	add_cat_comment(argc, argv, &thecat);
	write_cat_head(&thecat);
	while (read_object(&obj))
		write_object(&obj);
	fclose(opf);
	fclose(catf);
	if (0 > system("/bin/rm getpsf.tmp2\0"))
		error_exit("getpsf: system call 3 failed\n");

	/* now we analyse */
	if (makescript) {
		sprintf(systemstring, "analyse -f getpsf.tmp");
		fprintf(stderr, "# systemstring = %s\n", systemstring);
		if (0 > system(systemstring))
			error_exit("getpsf: system call 4 failed\n");
		if (!(catf = fopen("getpsf.tmp", "r")))
			error_exit("psfcorrect: failed to open getpsf.tmp\n");

		/* read the stars */
		set_cat_ipf(catf);
		read_cat_head(&thecat);
		nstars = 0;
		while (read_object(&obj)) {
			i[nstars] = obj.i;
			j[nstars] = obj.j;
			e1[nstars] = obj.e[0];
			e2[nstars] = obj.e[1];
			rh[nstars] = obj.rh;
			l[nstars] = (obj.l > 0 ? log(obj.l) : 0.0);
			nstars++;
			if (nstars >= MAX_STARS)
				error_exit("getpsf: too many stars selected; MAX_STARS = 1000\n");
		}
		fclose(catf);
		fprintf(stderr, "# %d stars read\n", nstars);
		if (0 > system("rm getpsf.tmp"))
			error_exit("getpsf: system call 5 failed\n");
	
		/* now make sub-sample in a box on e1, e2 plane */
		fprintf(stderr, "# getpsf: now you select a subset of the stars\n");
		fprintf(stderr, "# getpsf: same cursor rules as before\n");
		smmode = EEMODE;
		smpopup(drawfn, cursorfn, "");
		e1min = MIN(smx1, smx2);
		e1max = MAX(smx1, smx2);
		e2min = MIN(smy1, smy2);
		e2max = MAX(smy1, smy2);
		fprintf(stderr, "# selected range\n#\te1 = %10.3e to %10.3e\n#\te2 = %10.3e to %10.3e\n",
			e1min, e1max, e2min, e2max);
	
		/* take average ellipticity and rh of surviving subsample */
		ngood = 0;
		e1av = e2av = rav = 0.0;
		fprintf(stderr, "#    i     j         e1         e2          r\n");
		for (istar = 0; istar < nstars; istar++) {
			if (e1[istar] > e1max || e1[istar] < e1min || e2[istar] > e2max || e2[istar] < e2min)
				continue;
			fprintf(stderr, "#%5d %5d %10.3e %10.3e %10.3e\n", 
				i[istar], j[istar], e1[istar], e2[istar], rh[istar]);
			rav += rh[istar];
			e1av += e1[istar];
			e2av += e2[istar];
			ngood ++;
		}
		if (ngood) {
			e1av /= ngood;
       			e2av /= ngood;
			rav /= ngood;
			fprintf(stderr, "#\n# %d stars used; average properties:\n", ngood);
			fprintf(stderr, "#        e1         e2          r\n");
			fprintf(stderr, "#%10.3e %10.3e %10.3e\n", e1av, e2av, rav);
			e = sqrt(e1av * e1av + e2av * e2av);
			b = sqrt(rav * rav * (1 - e));
			a = sqrt(rav * rav / (1 - e));
			phi = 0.5 * atan2(e2av , e1av);
			fprintf(stderr, "#         a          b        phi\n");
			fprintf(stderr, "#%10.3e %10.3e %10.3e\n", a, b, 360 * phi / (2 * PI));
		} else {
			error_exit("getpsf: you fool, you selected no stars!\n");
		}
		fprintf(paramf, "%10.3e %10.3e %10.3e %s\n", a, b, 360 * phi / (2 * PI), catfilename);
		}
	}
	exit(0);
}


void	drawfn(void)
{
	float	dottype = 11.0, crosstype = 41.0, startype = 62.0;
	switch (smmode) {
		case RLMODE:
			sm_limits(0.0, 5.0, 3.0, 13.0);
			sm_box(1, 2, 0, 0);
			sm_ptype(&dottype, 1);
			sm_points(smr, sml, smN);
			sm_xlabel("r [pixels]");
			sm_ylabel("log(l)");
			break;
		case EEMODE:
			sm_limits(-ELIMIT, ELIMIT, -ELIMIT, ELIMIT);
			sm_box(1, 2, 0, 0);
			sm_ptype(&startype, 1);
			sm_points(e1, e2, nstars);
			sm_xlabel("e1");
			sm_ylabel("e2");
			sm_relocate(0.0, -ELIMIT);
			sm_draw(0.0, ELIMIT);
			sm_relocate(-ELIMIT, 0.0);
			sm_draw(ELIMIT, 0.0);
			break;
		default:
			error_exit("getpsf: bad smmode\n");
			break;
	}
}




void	cursorfn(float x, float y, int key)
{
	static	int first = 1;

	if (key == 'b') {
		if (first) {
			smx1 = x;
			smy1 = y;
			first = 0;
		} else {
			smx2 = x;
			smy2 = y;
			sm_erase();
			drawfn();
			drawselectionbox();
			first = 1;
		}
	}
}





void	drawselectionbox(void)
{
	sm_relocate((float) smx1, (float) smy1);
	sm_draw((float) smx1, (float) smy2);
	sm_draw((float) smx2, (float) smy2);
	sm_draw((float) smx2, (float) smy1);
	sm_draw((float) smx1, (float) smy1);
	sm_gflush();
}




