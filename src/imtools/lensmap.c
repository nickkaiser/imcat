/*
 * lensmap.c
 */

#define	usage "\n\n\n\
NAME\n\
	lensmap --- map a source image to a target image\n\
\n\
SYNOPSIS\n\
	lensmap	[option...]\n\
		-d  datafile	# file for cluster parameters (cluster.dat)\n\
		-o  xo, yo		# offset for the origin of the image (0,0)\n\
		-m mode				# mapping mode (1)\n\
\n\
DESCRIPTION\n\
	\"lensmap\" maps a source image to a target image using multiple\n\
	isothermal sphere lens model\n\
	lens properties specified in datafile in form\n\
\n\
	rEinstein rCore xc yc\n\
\n\
	all in pixel units\n\
\n\
	Use -m option to specify mode, where these are (in order of expense)\n\
		mode = 0:	# nearest pixel\n\
		mode = 1:	# linear interpolation\n\
		mode = 2:	# sum over triangles\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"		
	

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../utils/error.h"
#include "../imlib/fits.h"
#include "../imlib/map.h"
#include "lensmap.h"

#define  MAX_CLUSTERS	1000

static	float		re[MAX_CLUSTERS], rc[MAX_CLUSTERS];
static	int		nclusters, ic[MAX_CLUSTERS], jc[MAX_CLUSTERS], io, jo;

main(int argc, char *argv[])	
{
	int		arg = 1; 
	FILE		*clusterf;
	fitsheader	*fits;
	char		clusterfilename[128], inputline[1024];
	int		nc, N1, N2, mapmode, pixtype;
	float		**fsource, **ftarget;
	
	/* defaults */
	io = jo = 0;
	strcpy(clusterfilename, "cluster.dat");
	mapmode = FAST_MAP_MODE;

	while (arg < argc) {
		if (*argv[arg] == '-') {
			switch (*(argv[arg++]+1)) {
				case 'd':
					if (1 != sscanf(argv[arg++], "%s", clusterfilename))
						error_exit(usage);
					break;
				case 'o':
					if (1 != sscanf(argv[arg++], "%d", &jo))
					if (1 != sscanf(argv[arg++], "%d", &jo))
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%d", &io))
						error_exit(usage);
					break;
				case 'm':
					if (1 != sscanf(argv[arg++], "%d", &mapmode))
						error_exit(usage);
					break;
				default:
					error_exit(usage);
					break;
			}
		} else {
			error_exit(usage);
		}
	}
	
	/* read the cluster data */
	if (!(clusterf = fopen(clusterfilename, "r")))
		error_exit("simulatecat: failed to open cluster data file\n");
	nc = 0;
	while (fgets(inputline, 1024, clusterf)) {
		if (inputline[0] == '#')
			continue;
		if (4 != sscanf(inputline, "%f %f %d %d", re + nc, rc + nc, jc + nc, ic + nc))
			error_exit("simulatecat: bad format cluster data file\n");
		nc++;
	}
	nclusters = nc;
	for (nc = 0; nc < nclusters; nc++) {
		fprintf(stderr, "cluster props: re = %10.3f; rc = %10.3f; yc = %4d; xc = %4d\n",
			re[nc], rc[nc], ic[nc], jc[nc]);
		ic[nc] -= io;
		jc[nc] -= jo;
	}

	/* read the source file */
	read2Dfloatimage(&fsource, &N1, &N2, &fits, stdin);
	
	/* create the target image */
	allocFloatArray(&ftarget, N1, N2);
	
	/* apply the mapping */
	switch(mapmode) {
		case ULTRAFAST_MAP_MODE:
			ultrafastmap(ftarget, N1, N2, fsource, N1, N2, deflection);
			break;
		case FAST_MAP_MODE:
			fastmap(ftarget, N1, N2, fsource, N1, N2, deflection);
			break;
		case TRIANGLE_MAP_MODE:
			map(ftarget, N1, N2, fsource, N1, N2, deflection);
			break;
		default:
			error_exit("transformimage: bad mapping mode\n");
			break;
	}
	
	/* write the target image */
	add_comment(argc, argv, fits);
	write2Dfloatimage(ftarget, fits);
	
	/* all done */
	exit(0);
}



int	deflection(float ri, float rj, float *di, float *dj)
{
	float	si, sj, s, ss;
	int	nc;
	
	*di = *dj = 0;
	for (nc = 0; nc < nclusters; nc++) {
		si = ri - ic[nc];
		sj = rj - jc[nc];
		ss = si * si + sj * sj;
		if (ss > 0.0) {
			s = sqrt(ss);
			*di -= re[nc] * si / (s * (1 + (rc[nc] * rc[nc]) / ss));
			*dj -= re[nc] * sj / (s * (1 + (rc[nc] * rc[nc]) / ss));
		}
	}
	return(1);
}




