/*
 * readusnoacat.c
 */

#define usage "\n\
NAME\n\
	readusnoacat --- extract catalogue from USNO-A database\n\
\n\
SYNOPSIS\n\
	readusnoacat ra dl dec ddec\n\
\n\
DESCRIPTION\n\
	Readusnoacat extracts a lc-format catalogue from the\n\
	US Naval Observatory all-sky astrometric catalogue.\n\
\n\
	All objects with ra lying in the range\n\
	ra +- (dl / cos(dec)) and dec in the range\n\
	dec +- ddec (epoch 2000) are extracted. \n\
	The parameterisation of the range in longitude dl\n\
	is chosen so that in the small angle approximation\n\
	one obtains a box of height ddec and width dl both in\n\
	degrees.\n\
\n\
	Angle arguments may be given in decimal notation, in\n\
	which case they are interpreted as degrees, or as\n\
	colon separated triplets, in which case they are interpreted\n\
	as h:m:s (for ra, dra) and d:m:s (dec, ddec)\n\
\n\
	Readusnoacat expects to find an environment variable\n\
	USNOADIR telling it the directory containing the source\n\
	catalogue files.\n\
\n\
	The output catalogue contains the following entries:\n\
		x[2]		# sterographic sky coords\n\
		RA		# right ascension [deg]\n\
		DEC		# declination [deg]\n\
		rmag		# red magnitude\n\
		bmag		# blue magnitude\n\
		qflag		# indicates mag error if set\n\
		gflag		# indicates a correlated GSC entry\n\
		fieldno		# source image number\n\
	See the USNO-A readme files for more information on the\n\
	meaning of the flags.\n\
\n\
	x[2] is a stereographic projection (in units of degrees).\n\
\n\
AUTHOR\n\
	Nick Kaiser	kaiser@hawaii.edu\n\
\n\n"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "utils/error.h"
#include "imlib/fits.h"
#include "readusnoacat.h"
#include "getxfromradecfunc.h"
#include "radecio.h"

#define	SPDBANDWIDTH10 75
#define	SPDBANDWIDTH 7.5

#define swap_word(a) ( ((a) << 24) | \
                      (((a) << 8) & 0x00ff0000) | \
                      (((a) >> 8) & 0x0000ff00) | \
        ((unsigned int)(a) >>24) )

FILE	*lcpipe;

main(int argc, char *argv[])
{
	int	count, arg;
	double	ra, dec, dl, ddec;
	double	ra0, ra1, dec0, dec1, spd0, spd1, spdmin, spdmax;
	char	*usnaodatadir, basefilename[1024];
	int	ispd;
	char	lcstring[64], argstring[512];
	int	type;

	/* parse args */
	if (argc < 5) {
		error_exit(usage);
	}
	arg = 1;
	ra = getra(argv[arg++], &type);		
	dl = getra(argv[arg++], &type);
	dec = getdec(argv[arg++], &type);
	ddec = getdec(argv[arg++], &type);
	dec0 = dec - ddec;
	dec1 = dec + ddec;
	if ((dec0 < -90.0) || (dec1 > 90.0)) {
		error_exit("readusnoacat: illegal dec range - bailing\n");
	}		
	dl /= cos(M_PI * dec / 180.0);
	ra0 = ra - dl;
	ra1 = ra + dl;
	if ((ra0 < 0.0) || (ra1 > 360.0)) {
		error_exit("readusnoacat: illegal ra range - bailing\n");
	}
		

	/* get the data directory env variable */
	usnaodatadir = getenv("USNOADIR");
	if (!usnaodatadir) {
		error_exit("readusnoacat: can't source catalogue; please define env variable USNOADIR\n");
	}

	/* open the output pipe */
	argsToString(argc, argv, argstring);
	sprintf(lcstring, "lc -C -x -a \"history: %s\" -N '1 2 x' -n RA -n DEC -n rmag -n bmag -n qflag -n gflag -n fieldno", argstring);
	lcpipe = popen(lcstring, "w");
        if (!lcpipe) {
                fprintf(stderr, "readusnoacat: unable to open lcpipe for output\n");
                exit(-1);
        }

	/* figure out which source catalogues we need and extract objects*/
	count = 0;
	spdmax = dec1 + 90;
	spdmin = dec0 + 90; 
	for (ispd = 0; ispd < 24; ispd++) {
		spd0 = ispd * SPDBANDWIDTH;
		spd1 = (ispd + 1) * SPDBANDWIDTH;
		if ((spdmin < spd1) && (spdmax > spd0)) {
			sprintf(basefilename, "%s/zone%04d", usnaodatadir, ispd * SPDBANDWIDTH10);
			count +=getobjects(basefilename, ra0, ra1, dec0, dec1, ra, dec);
		} 
	}
	
	exit(pclose(lcpipe));
}


int 	getobjects(char *basefilename, double ra0, double ra1, double dec0, double dec1, double rao, double deco)
{
	FILE	*catfile, *accfile;
	char	catfilename[1024], accfilename[1024], accelrec[512];
	int	count = 0;
	double	ramin, ramax, ra, dec, x[2];
	long	offset;
	int	catfileisopen;
	int	rec[3], iRA, iDEC, iDAT;
	double	rmag, bmag;
	int	fieldno, qflag, gscflag;

	/* open the accelerator file */	
	sprintf(accfilename, "%s.acc", basefilename);
	accfile = fopen(accfilename, "r");
	if (!accfile) {
		fprintf(stderr, "readusnoacat: unable to open accelerator file %s\n", accfilename);
		exit(-1);
	}

	/* find segments we need */
	catfileisopen = 0;
	while(fgets(accelrec, 512, accfile)) {
		sscanf(accelrec, "%lf %ld", &ramin, &offset);
		/* offset is like a fortran index so decrement it */
		offset--;
		/* and convert to bytes */
		offset *= 12;
		ramax = ramin + 0.25;
		/* convert ramax, ramin to degrees */
		ramin *= 15.0;
		ramax *= 15.0;
		if ((ra0 < ramax) && (ra1 > ramin)) {
			/* open cat file and spool to position if necessary */
			if (!catfileisopen) {
				sprintf(catfilename, "%s.cat", basefilename);
				catfile = fopen(catfilename, "r");
				if (!catfile) {
					fprintf(stderr, "readusnoacat: unable to open catalogue file %s\n", catfilename);
					exit(-1);
				}
				fseek(catfile, offset, SEEK_SET);
			}
			while (fread(rec, 3, sizeof(int), catfile)) {
#ifdef LITTLE_ENDIAN
				iRA  = swap_word(rec[0]);
				iDEC = swap_word(rec[1]);
				iDAT = swap_word(rec[2]);
#else
				iRA  = rec[0];
				iDEC = rec[1];
				iDAT = rec[2];
#endif
				ra = iRA / (3600 * 100.0);
				dec = iDEC / (3600 * 100.0) - 90;
	
				if ((ra > ra0) && (ra < ra1) && (dec > dec0) && (dec < dec1)) {
					count++;
					/* decompose the data field */
					if (iDAT < 0) {
						gscflag = 1;
						iDAT *= -1;
					} else {
						gscflag = 0;
					}
					rmag = fmod((double) iDAT, 1000.0);
					iDAT -= (int) rmag;
					rmag /= 10.0;
					bmag = fmod((double) iDAT, 1000000.0);
					iDAT -= (int) bmag;
					bmag /= 10000.0;
					fieldno = (int) fmod((double) iDAT, 1000000000.0);
					iDAT -= fieldno;
					fieldno /= 1000000;
					qflag = (iDAT == 0 ? 0 : 1);
					getxcoords(ra, dec, rao, deco, x);
					fprintf(lcpipe, "%13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %d %d %d\n", x[0], x[1], ra, dec, rmag, bmag, qflag, gscflag, fieldno);
				}
		 	}
		}
	}

	close(catfile);
	close(accfile);

	return (count);
}


