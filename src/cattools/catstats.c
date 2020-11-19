/*
 * catstats.c
 */

#define usage "\n\n\
NAME\n\
	catstats - calculate statistics for catalogue object values\n\
\n\
SYNOPSIS\n\
	catstats [options....]\n\
		-s		# only output count, min, max, mean, sigma1\n\
		-v statistic	# only output 'statistic'\n\
		-o		# output one line per variable (tab spaced)\n\
\n\
DESCRIPTION\n\
	'catstats' reads a catalogue from stdin and creates a new\n\
	catalogue whose first object item is a text entry 'statistic'\n\
	which takes values 'count', 'min', 'max', 'mean'....\n\
	and whose subsequent object items have the same names as the\n\
	numerical items in the input catalogue and which\n\
	contain the appropriate statistic.\n\
\n\
	Catstats will always output the basic statistics:\n\
		count		# number of objects\n\
		min		# minimum value\n\
		max		# maximum value\n\
		mean		# <f> = sum f / count\n\
		sigma1		# sqrt(<f^2> - <f>^2)\t\n\
	and by default will also calculate\n\
		mode		# 'robust' mode estimator\n\
		median		# median\n\
		lquart		# upper quartile\n\
		uquart		# lower quartile\n\
		sigma2		# 'robust' sigma estimator   \n\
	provided there are are least 8 objects in the catalogue.\n\
\n\
	The statistics 'mode' and 'sigma2' are designed to be\n\
	insensitive to outliers.\n\
\n\
	The 'mode' is estimated by first making a crude\n\
	estimate of sigma as (uquart - lquart) / 1.34; smoothing\n\
	the histogram of values with a gaussian of width sigma\n\
	and returning the location of the peak.\n\
\n\
	'sigma2' is estimated from the width of the region around\n\
	the mode containing 25 percent of the values (and assuming)\n\
	a gaussian distribution.  Thus, sigma2 is effectively\n\
	measured from the curvature of the distribution around the\n\
	mode and will, for example, overestimate the real sigma\n\
	for a very boxy distribution.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n"


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "../utils/error.h"
#include "../utils/stats_stuff.h"
#include "../catlib/cat.h"



#define	SIMPLE_STATS_MODE	0
#define ALL_STATS_MODE		1
#define	N_STATS			9
#define MIN_INDEX		0
#define MAX_INDEX		1
#define MEAN_INDEX		2
#define SIGMA1_INDEX		3
#define MODE_INDEX		4
#define MEDIAN_INDEX		5
#define LQUART_INDEX		6
#define UQUART_INDEX		7
#define SIGMA2_INDEX		8


char	statname[N_STATS][64] = {
	"min",
	"max",
	"mean",
	"sigma1",
	"mode",
	"median",
	"lquart",
	"uquart",
	"sigma2"
};

typedef struct statobj {
	float		*f;
	struct statobj	*next;
} statobj;


void	getnumericaddresses(object *obj);
void	getitemaddresses(item *theitem, void *addr, int level);

#define MAX_NUMBERS	10000
#define MAX_VARNAMELEN	100
#define MIN_COUNT	8

static double	*fin[MAX_NUMBERS];
static double	*fstat[N_STATS][MAX_NUMBERS];
static int	g_istat, findex;
static char	g_varname[MAX_NUMBERS][MAX_VARNAMELEN];

main(int argc, char *argv[])
{
	int		arg = 1, opmode, nnumeric, i, iobj, singlestatistic, statindex;
	cathead		*ipcat, *opcat;
	object		*ipobj, *opobj[N_STATS];
	item		*ipitem, *opitem;
	int		istat, nstats, first, count;
	float		*flist, fmode, fmedian, fuquart, flquart, fsigma, *f[1], fjunk;
	statobj		*basestatobj = NULL, *newstatobj, *thestatobj;
	int		onelineoutput = 0;

	/* defaults */
	opmode = ALL_STATS_MODE;
	nstats = N_STATS;
	singlestatistic = 0;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch (argv[arg++][1]) {
			case 's':
				opmode = SIMPLE_STATS_MODE;
				nstats = 4;
				break;
			case 'o':
				onelineoutput = 1;
				break;
			case 'v':
				singlestatistic = 1;
				for (i = 0; i < N_STATS; i++) {
					if (!strcmp(argv[arg], statname[i])) {
						statindex = i;
						break;
					}
				}
				if (i == N_STATS) {
					error_exit("catstats: invalid argument following -v flag\n");
				}
				if (i < 4) {
					opmode = SIMPLE_STATS_MODE;
				}
				arg++;
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	/* read the cat header and input object*/
	ipcat  = readcathead();
	ipobj = newobject(ipcat);

	/* create the output cathead */
	opcat = (cathead *) calloc(1, sizeof(cathead));
	copyheaderinfo(opcat, ipcat);
	addargscomment(argc, argv, opcat);

	/* create the output object items */
	opitem = newitem("statistic", TEXT_TYPE, 1, 1);
	addobjectitem(opitem, opcat);
	ipitem = ipcat->objectitembase;
	i = 0;
	while (ipitem) {
		allocitemcontents(ipitem, &(ipitem->addr), 0);
		if (ipitem->itype == NUM_TYPE) {
			opitem = copyitem(ipitem);
			opitem->addr = ipitem->addr;
			addobjectitem(opitem, opcat);
		}
		ipitem = ipitem->next;
	}
	connectobjecttocathead(ipobj);

	/* write output cathead */
	if (!singlestatistic && !onelineoutput) {
		writecathead(opcat);
	}

	/* get all the numeric addresses of the input object as a single vector */
	g_istat = -1;
	getnumericaddresses(ipobj);
	nnumeric = findex;

	/* create the output objects */
	for (istat = 0; istat < nstats; istat++) {
		opobj[istat] = newobject(opcat);
		allocobjectcontents(opobj[istat]);
		((char **) ((opobj[istat])->addrlist)[0])[0] = statname[istat];
		g_istat = istat;
		getnumericaddresses(opobj[istat]);
	}

	first = 1;
	count = 0;
	while (readobject(ipobj)) {
		count++;
		for (i = 0; i < nnumeric; i++) {
			if (first) {
				*(fstat[MIN_INDEX][i]) = *(fin[i]);
				*(fstat[MAX_INDEX][i]) = *(fin[i]);
			} else {
				if (*(fin[i]) < *(fstat[MIN_INDEX][i]))
					*(fstat[MIN_INDEX][i]) = *(fin[i]);
				if (*(fin[i]) > *(fstat[MAX_INDEX][i]))
					*(fstat[MAX_INDEX][i]) = *(fin[i]);
			}
			*(fstat[MEAN_INDEX][i]) += *(fin[i]);
			*(fstat[SIGMA1_INDEX][i]) += *(fin[i]) * *(fin[i]);
		}
		first = 0;
		if (opmode == ALL_STATS_MODE) {
			newstatobj = (statobj *) calloc(1, sizeof(statobj));
			newstatobj->f = (float *) calloc(nnumeric, sizeof(float));
			for (i = 0; i < nnumeric; i++) {
				(newstatobj->f)[i] = (float) *(fin[i]);
			}
			newstatobj->next = basestatobj;
			basestatobj = newstatobj;
		}
	}
	flist = (float *) calloc(count, sizeof(float));
     
	for (i = 0; i < nnumeric; i++) {
		*(fstat[MEAN_INDEX][i]) /= count;
		*(fstat[SIGMA1_INDEX][i]) /= count;
		*(fstat[SIGMA1_INDEX][i]) -= *(fstat[MEAN_INDEX][i]) * *(fstat[MEAN_INDEX][i]);
		if (*(fstat[SIGMA1_INDEX][i]) > 0.0)
			*(fstat[SIGMA1_INDEX][i]) = sqrt(*(fstat[SIGMA1_INDEX][i]));
		else
			*(fstat[SIGMA1_INDEX][i]) = 0.0;
		if (opmode == ALL_STATS_MODE && count >= MIN_COUNT) {
			iobj = 0;
			thestatobj = basestatobj;
			while (thestatobj) {
				flist[iobj++] = (thestatobj->f)[i];
				thestatobj = thestatobj->next;
			}
			liststats(flist, count, &fmedian, &flquart, &fuquart, &fsigma);
			f[0] = flist;
			findmode(f, count, 1, flquart, fuquart, &fmode, &fjunk);	
			*(fstat[MODE_INDEX][i]) = (double) fmode;
			*(fstat[MEDIAN_INDEX][i]) = (double) fmedian;
			*(fstat[LQUART_INDEX][i]) = (double) flquart;
			*(fstat[UQUART_INDEX][i]) = (double) fuquart;
			*(fstat[SIGMA2_INDEX][i]) = (double) fsigma;
		}
	}

	if (singlestatistic) {
		for (i = 0; i < nnumeric; i++) {
			fprintf(stdout, "%.8lg", *(fstat[statindex][i]));
			if (i + 1 < nnumeric) {
				fprintf(stdout, " ");
			}
		}
		fprintf(stdout, "\n");
		exit(0);
	}

	if (onelineoutput) {
		for (i = 0; i < nnumeric; i++) {
			if (opmode == SIMPLE_STATS_MODE) {
				nstats = N_STATS;
			}
			fprintf(stdout, "%s\t", g_varname[i]);
			for (statindex = 0; statindex < nstats; statindex++) {
				fprintf(stdout, "%13.8lg", *(fstat[statindex][i]));
				if (statindex < (nstats - 1)) {
					fprintf(stdout, "\t");
				}
			}
			fprintf(stdout, "\n");
		}
	} else {
		for (istat = 0; istat < nstats; istat++) {
			writeobject(opobj[istat]);
		}
	}
	exit(0);	
}



void	getnumericaddresses(object *obj)
{
	item	*theitem;

	findex = 0;
	theitem = (obj->cathead)->objectitembase;
	while (theitem) {
		if (theitem->itype == NUM_TYPE) {
			getitemaddresses(theitem, theitem->addr, 0);
		}
		theitem = theitem->next;
	}
}

void	getitemaddresses(item *theitem, void *addr, int level)
{
	int	i;

	if (level < theitem->ndim - 1) {
		for (i = 0; i < (theitem->dim)[level]; i++) {
			getitemaddresses(theitem, *((void **) addr + i), level + 1);
		}
	} else {
		for (i = 0; i < (theitem->dim)[level]; i++) {
			if ((theitem->dim)[level] > 1) {
				snprintf(g_varname[findex], MAX_VARNAMELEN, "%s[%d]", theitem->name, i);
			} else {
				snprintf(g_varname[findex], MAX_VARNAMELEN, "%s", theitem->name);
			}
			if (g_istat < 0)
				fin[findex++] = (double *) addr + i;
			else
				fstat[g_istat][findex++] = (double *) addr + i;
		}
	}
}

