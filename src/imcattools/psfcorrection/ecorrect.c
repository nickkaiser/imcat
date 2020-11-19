#define	usage "\n\n\n\
NAME\n\
	ecorrect --- correct ellipticities for psf anisotropy\n\
\n\
SYNOPSIS\n\
	ecorrect [option...]\n\
		-f  efitdata		# file for stellar e / psm model parameters (efit.out)\n\
\n\
DESCRIPTION\n\
	\"ecorrect\" reads a catalogue from stdin and corrects ellipticities\n\
	according to model for psf polarization p = e / psm as determined by efit.\n\
	Correction applied is e -= psm * p.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"		
	

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "../../utils/error.h"
#include "../../catlib/cat.h"
#include "efit_stuff.h"

main(int argc, char *argv[])	
{
	int		arg = 1, shift, i; 
	FILE		*efitdataf;
	int		comc, mode, nmodes, N, order;
	char		*efitfilename, defefitfilename[10] = "efit.out", line[2048];
	double		**amp, p0, p1, *x, *e, **psm;
	cathead		*thecathead;
	object		*theobject;
	
	
	
	/* defaults */
	efitfilename = defefitfilename;
	
	while (arg < argc) {
		if (*argv[arg] == '-') {
			switch (*(argv[arg++]+1)) {
				case 'f':
					efitfilename = argv[arg++];
					break;
				default:
					error_exit(usage);
					break;
			}
		} else {
			error_exit(usage);
		}
	}
	
	/* read the efit data */
	if (!(efitdataf = fopen(efitfilename, "r")))
		error_exit("ecorrect: failed to open e-fit data file\n");
	readamplitudes(&nmodes, &N, &order, &amp, efitdataf);
	fprintf(stderr, " %d mode model\n# %d order\n# %d pixels wide\n",
		nmodes, order, N);
	setframesize(N);
	
	setcatopfiletype(BINARY_FILE_TYPE);

	thecathead = readcathead();
	theobject = newobject(thecathead);
	connectobjecttocathead(theobject);
	allocobjectcontents(theobject);
	x = (double *) ((theobject->addrlist)[getobjectitemindex("x", theobject)]);
	e = (double *) ((theobject->addrlist)[getobjectitemindex("e", theobject)]);
	psm = (double **) ((theobject->addrlist)[getobjectitemindex("psm", theobject)]);
	addargscomment(argc, argv, thecathead);
	writecathead(thecathead);
	while (readobject(theobject)) {
		p0 = p1 = 0;
		for (mode = 0; mode < nmodes; mode++) {
			p0 += amp[0][mode] * g(mode, x[0], x[1]);
			p1 += amp[1][mode] * g(mode, x[0], x[1]);
		}
		e[0] -= (psm[0][0] * p0 + psm[0][1] * p1);
		e[1] -= (psm[1][0] * p0 + psm[1][1] * p1);
		writeobject(theobject);
	}
	exit(0);
}



