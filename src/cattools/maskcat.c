/*
 *	maskcat.c 
 */

#define	usage "\n\n\n\
NAME\n\
	maskcat -- remove objects falling in a list of rectangles\n\
SYNOPSIS\n\
	maskcat	[options....] maskfile \n\
		-x xname	# name for spatial coordinate ('x')\n\
		-i		# inverse selection mode\n\
\n\
DESCRIPTION\n\
	\"maskcat\" reads a catalogue from stdin, removes any objects\n\
	which fall within any of a set of rectangles specified in the file\n\
	'maskfile', and writes the result to stdout.\n\
\n\
	The mask file must be a lc-format catalogue containing\n\
	(at least) a pair of position vectors x1[2], x2[2]\n\
	which are diagonally opposite corners of the rectangle.\n\
	Objects lying on the boundary will be rejected.\n\
\n\
	With the -i option, we output only those objects which fall\n\
	within the union of the rectangles sepcified in maskfile.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"		


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../catlib/cat.h"					/* READ THIS! */
#include "../utils/error.h"


#define MAX_RECTS	10000


main(int argc, char *argv[])	
{
	int		arg = 1, nrects, irect, inrect;
	double		x1[MAX_RECTS], y1[MAX_RECTS], x2[MAX_RECTS], y2[MAX_RECTS], temp;
	char		tempstring[1024], *xname, defxname[128] = "x", line[1024];
	FILE		*maskf;
	cathead		*inputcathead, *outputcathead;
	object		*theobject;
	double		*x;
	int		inverseselection;

	/* defaults */
	xname = defxname;
	inverseselection = 0;
	
	while (arg < (argc - 1)) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		} else {
			switch (argv[arg++][1]) {
				case 'x':
					xname = argv[arg++];
					break;
				case 'i':
					inverseselection = 1;
					break;
				default:
					error_exit(usage);
					break;
			}
		}
	}
	if (argv[arg][0] == '-')
		error_exit(usage);

	sprintf(tempstring, "lc -o x1 x2 < %s", argv[arg]);

	maskf = popen(tempstring, "r");
	if (!maskf)
		error_exit("maskcat: failed to open mask file\n");
	nrects = 0;
	while (fgets(line, 1024, maskf)) {
		if (line[0] == '#')
			continue;
		if (4 != sscanf(line, "%lf %lf %lf %lf", x1 + nrects, y1 + nrects, x2 + nrects, y2 + nrects)) {
			error_exit("maskcat: corrupted mask file?\n");
		}
		if (x1[nrects] > x2[nrects]) {
			temp = x1[nrects];
			x1[nrects] = x2[nrects];
			x2[nrects] = temp;
		}
		if (y1[nrects] > y2[nrects]) {
			temp = y1[nrects];
			y1[nrects] = y2[nrects];
			y2[nrects] = temp;
		}
		nrects++;
		if (nrects == MAX_RECTS)
			error_exit("maskcat: too many rectangles\n");
	}
	pclose(maskf);

	setcatopfiletype(BINARY_FILE_TYPE);

	inputcathead = readcathead();			/* read the cat head */
	theobject = newobject(inputcathead);		/* make the input object */
	connectobjecttocathead(theobject);		/* obj addresses point back to cathead */
	allocobjectcontents(theobject);		/* and allocate space for obj data */

	outputcathead = (cathead *) calloc(1, sizeof(cathead));			/* new cathead */
	copyheaderinfo(outputcathead, inputcathead);				/* copy header stuff */

	addargscomment(argc, argv, outputcathead);		/* add history */

	copycontentinfo(outputcathead, inputcathead);		/* copy over pre-exisiting object items */
	writecathead(outputcathead);				/* and write cathead out */			

	x = (double *) ((theobject->addrlist)[getobjectitemindex(xname, theobject)]);

	while (readobject(theobject)) {					/* big loop */
		inrect = 0;
		for (irect = 0; irect < nrects; irect++) {
			if (x[0] >= x1[irect] && x[0] <= x2[irect] && x[1] >= y1[irect] && x[1] <= y2[irect]) {
				inrect = 1;
				break;
			}
		}
		if (inverseselection) {
			if (inrect) {
				writeobject(theobject);
			}
		} else {
			if (!inrect) {
				writeobject(theobject);
			}
		}
	}
	exit(0);
}

