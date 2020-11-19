#define	usage "\n\n\n\
SYNOPSIS\n\
		ecorrect	[option...] a.cat b.cat .....\n\
			-f  efitdata		# file for stellar e / psm model parameters (efit.out)\n\
\n\
DESCRIPTION\n\
		\"ecorrect\" corrects ellipticities in number of catalogues\n\
		according to model for psf polarization p = e / psm as determined by efit\n\
		correction applied is e -= psm * p. If you use multiple cats\n\
		be sure that they correspond precisely to those fed to efit.\n\
		Ecorrect generates a new catalogue with suffix '.cat.x' for\n\
		each input cat.\n\
\n\n\n"		
	

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "../../utils/error.h"
#include "../../catlib/cat.h"
#include "efit_stuff.h"

#define MAX_CATS 1000

main(int argc, char *argv[])	
{
	int		arg = 1, shift, i; 
	FILE		*efitdataf, *catf, *opf;
	int		comc, mode, nmodes, frame, nframes, ncats = 0, N, order;
	char		efitdatafilename[128], systemstring[1024], catfilename[MAX_CATS][128];
	char		tempfile[128], format[128] = "", line[2048];
	double		**amp, p0, p1, *x, *e, **psm;
	cathead		*thecathead;
	object		*theobject;
	
	
	
	/* defaults */
	strcpy(efitdatafilename, "efit.out");
	
	while (arg < argc) {
		if (*argv[arg] == '-') {
			switch (*(argv[arg++]+1)) {
				case 'f':
					if (EOF == sscanf(argv[arg++], "%s", efitdatafilename))
						error_exit(usage);
					break;
				default:
					error_exit(usage);
					break;
			}
		} else {
			strcpy(catfilename[ncats++], argv[arg++]);
			if (ncats >= MAX_CATS)
				error_exit("ecorrect: too many cat files\n");
		}
	}
	
	/* read the efit data */
	if (!(efitdataf = fopen(efitdatafilename, "r")))
		error_exit("ecorrect: failed to open e-fit data file\n");
	readamplitudes(&nmodes, &nframes, &N, &order, &amp, efitdataf);
	if (nframes != ncats)
		error_exit("ecorrect: mismatch between nframes and ncats given as args\n");
	fprintf(stderr, " %d mode model\n# %d order\n# %d frames\n# %d pixels wide\n",
		nmodes, order, nframes, N);
	setnframes(nframes);
	setframesize(N);
	
	sprintf(tempfile, "%d.tmp", getpid());
	for (i = 0; i < 8; i++) {
		strcat(format, LC_NUM_FMT);
	}

	setcatopfiletype(BINARY_FILE_TYPE);

	for (frame = 0; frame < nframes; frame++) {
		if (!(catf = fopen(catfilename[frame], "r")))
			error_exit("ecorrect: failed to open lc pipe for input\n");		
		if (!(opf = fopen(tempfile, "w")))
			error_exit("ecorrect: failed to open temp-file for output\n");
		setcatipf(catf);
		setcatopf(opf);
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
				p0 += amp[0][mode] * g(mode, frame, x[0], x[1]);
				p1 += amp[1][mode] * g(mode, frame, x[0], x[1]);
			}
			e[0] -= (psm[0][0] * p0 + psm[0][1] * p1);
			e[1] -= (psm[1][0] * p0 + psm[1][1] * p1);
			writeobject(theobject);
		}
		sprintf(systemstring, "mv %s %s.x", tempfile, catfilename[frame]);
		system(systemstring);
	}
	exit(0);
}



