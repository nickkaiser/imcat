#define	usage "\n\n\n\
NAME\n\
	hcat2cat - process hfindpeaks output\n\
\n\
SYNOPSIS\n\
	hcat2cat [options...] < hcatfile > catfile\n\
		-a		# fire on all local maxima\n\
DESCRIPTION\n\
	\"hcat2cat\" reads a set of peak trajectories in \"hcat\"\n\
	format from stdin and applies an algorithm to pick out\n\
	particular points (e.g. points of max significance).\n\
	Default is to pick only the most significant local maxima\n\
	along a peak trajectory, but -a option finds all local maxima.\n\
	Standard format catalogue goes to stdout\n\
	THIS IS NOW CALLED AUTOMATICALLY BY hfindpeaks\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"		

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../../utils/error.h"
#include "../../utils/stats_stuff.h"
#include "../../catlib/cat.h"
#include "hfindpeaks.h"

typedef struct oldobject {
        int             i, j;                   /* peak position */
        float   nu, fs;                         /* peak parameters */
        float   rh, rg;				/* radius assigned by hcat2cat */
        float   l;                              /* luminosity */
        float   lg;                             /* hfindpeaks luminosity */
        float   e[2];                           /* ellipticity components */
        float	x[2];				/* accurate position */
} oldobject;


#define BESTPEAK	0
#define ALLPEAKS	1

void	process(peak *thepeak, int peakcount, int mode);
float	interp(float d, float x1, float x2, float x3);
void	write_object(oldobject *objectptr);

/* empirical factors to convert rg => rh, fs => mag */
#define RFACTOR         0.66
#define FFACTOR         15.41
#define EFACTOR         1.50

static object	*theobject;
static double	X[2], LG, RG, EG[2], FS, NU;

main(int argc, char *argv[])	
{
	cathead		*thecathead;
	line		*linehead = NULL, *theline;
	peak		*newpeak, *thepeak;
	int		arg = 1, peakcount = 0, mode;

	/* defaults */
	mode = BESTPEAK;

	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch (argv[arg++][1]) {
			case 'a':
				mode = ALLPEAKS;
				break;
			default:
				error_exit(usage);
		}
	}
	
	setcatopfiletype(BINARY_FILE_TYPE);
	thecathead = readcathead();
	addargscomment(argc, argv, thecathead);
	theobject = newobject(thecathead);
	connectcatheadtoobject(theobject);
        setaddress(theobject, getobjectitemindex("x", theobject), X);
        setaddress(theobject, getobjectitemindex("fs", theobject), &FS);
        setaddress(theobject, getobjectitemindex("nu", theobject), &NU);
        setaddress(theobject, getobjectitemindex("eg", theobject), EG);
        setaddress(theobject, getobjectitemindex("lg", theobject), &LG);
        setaddress(theobject, getobjectitemindex("rg", theobject), &RG);
	writecathead(thecathead);
	theline = linehead;
	newpeak = (peak *) calloc(1, sizeof(peak));
	thepeak = NULL;
	while(1 == fread(newpeak, sizeof(peak), 1, stdin)) {
		peakcount++;
		if (newpeak->next == NULL) {	/* end of a line */
			newpeak->next = thepeak;
			thepeak = newpeak;
			process(thepeak, peakcount, mode);
			peakcount = 0;
			thepeak = NULL;
		} else {
			newpeak->next = thepeak;
			thepeak = newpeak;
		}
		newpeak = (peak *) calloc(1, sizeof(peak));
	}
	exit(0);
}




void	process(peak *thepeak, int peakcount, int mode)
{
	peak		*nextpeak, **pk;
	int		i = 0, count = 0;
	static	oldobject	*obj = NULL, *bestobj, *tempobj;
	float		maxnu = 0.0, d, y1, y2, y3, bestnu = 0.0;

	if (!obj)
		obj = (oldobject *) calloc(1, sizeof(oldobject));
	if (!bestobj)
		bestobj = (oldobject *) calloc(1, sizeof(oldobject));

	pk = (peak **) calloc(peakcount, sizeof(peak));

	while (thepeak) {
		pk[i++] = thepeak;
		thepeak = thepeak->next;
	}

	for (i = 1; i < peakcount - 1; i++) {
		if (pk[i]->nu > pk[i-1]->nu && pk[i]->nu > pk[i+1]->nu) {
			y1 = pk[i-1]->nu;
			y2 = pk[i]->nu;
			y3 = pk[i+1]->nu;
			d = 0.5 * (y3 - y1) / (2 * y2 - y1 - y3);
			d = (d > 1 ? 1 : d);
			d = (d < -1 ? -1 : d);
			obj->i = pk[i]->i;
			obj->j = pk[i]->j;
			obj->x[0] = pk[i]->x[0];
			obj->x[1] = pk[i]->x[1];
			obj->e[0] = EFACTOR * pk[i]->e1;
			obj->e[1] = EFACTOR * pk[i]->e2;
			obj->rg = interp(d, pk[i-1]->rf, pk[i]->rf, pk[i+1]->rf);
			obj->fs = interp(d, pk[i-1]->fs, pk[i]->fs, pk[i+1]->fs);
			obj->nu = interp(d, pk[i-1]->nu, pk[i]->nu, pk[i+1]->nu);
			obj->l = obj->lg = FFACTOR * pk[i]->fs * pk[i]->rf * pk[i]->rf;
			obj->rh = RFACTOR * obj->rg;
			if (mode == ALLPEAKS) {
				write_object(obj);
				continue;
			} else {
				if (obj->nu > bestnu) {
					bestnu = obj->nu;
					/* swap them */
					tempobj = obj;
					obj = bestobj;
					bestobj = tempobj;
				}
			}
		}
	}

	if (mode == BESTPEAK && bestnu > 0) 
		write_object(bestobj);

	for (i = 0; i < peakcount; i++)
		free(pk[i]);
	free(pk);
}


float	interp(float d, float x1, float x2, float x3)
{
	if (d > 0)
		return((1 - d) * x2 + d * x3);
	else
		return((1 + d) * x2 - d * x1);
}



void	write_object(oldobject *obj)
{
	X[0] = (double) obj->x[0];
	X[1] = (double) obj->x[1]; 
	LG = (double) obj->lg;
	RG = (double) obj->rg;
	EG[0] = -(double) obj->e[0];
	EG[1] = (double) obj->e[1]; 
	FS = (double) obj->fs;
	NU = (double) obj->nu;
	writeobject(theobject);
}


