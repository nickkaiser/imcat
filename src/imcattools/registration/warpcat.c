/*
 * warpcat.c
 */

#define	usage "\n\n\n\
NAME\n\
	warpcat - apply spatial transformation to a catalogue\n\
\n\
SYNOPSIS\n\
	warpcat [options...] distparfile\n\
		-p psi_xx psi_xy psi_yx psi_yy\n\
		-d d_x d_y\n\
\n\
DESCRIPTION\n\
	\"warpcat\" reads a catalogue from stdin and applies a\n\
	spatial transformation to position vector x[2] according\n\
	to the parameters in 'distparfile'. It sends to stdout\n\
	a catalogue with extra vector\n\
		r = x + sum_m a_m f_m(x)\n\
	where mode functions are polynomials in x[0], x[1].\n\
	Use -p and -d options to apply a further linear transformation:\n\
		r_i => phi_ij r_j + d_j\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@ifa.hawaii.edu\n\
\n\n\n"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "utils/error.h"
#include "catlib/cat.h"
#include "imlib/fits.h"
#include "utils/modefunc.h"

void		rmodel(double *x, double *r);
float		phi[2][2], d[2];
int		nmodes, *l, *m;
double		**a;

#define TINY	1.e-50

main(int argc, char *argv[])	
{
	int		comc, arg = 1;
	char		line[1024];
	FILE		*distparfile;
	int		mode, cattype;
        cathead         *inputcathead, *outputcathead;
        object          *inputobject, *outputobject;
	item		*ritem;
	int		rindex;
	int		i, dodisplacement, dophitrans, dodistortions = 1, ndim;
	double		*x, *r, rp[2];
	char		*vardef[MODEFUNC_MAX_VARS], *xname;
	int		nvar;

	/* defaults */
	dodisplacement = 0;
	dophitrans = 0;

	/* parse args */
	if (argc < 2)
		error_exit(usage);
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			break;
		} else {
			switch (argv[arg++][1]) {
				case 'd':
					dodisplacement = 1;
					if ((argc - arg) < 2)
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%f", &(d[0])) ||
						1 != sscanf(argv[arg++], "%f", &(d[1])))
						error_exit(usage);
					break;
				case 'p':
					dophitrans = 1;
					if ((argc - arg) < 4)
						error_exit(usage);
					if (1 != sscanf(argv[arg++], "%f", &(phi[0][0])) ||
						1 != sscanf(argv[arg++], "%f", &(phi[0][1])) ||
						1 != sscanf(argv[arg++], "%f", &(phi[1][0])) ||
						1 != sscanf(argv[arg++], "%f", &(phi[1][1])))
						error_exit(usage);
					break;
				case 'u':
				default:
					error_exit(usage);
					break;
			}
		}
	}

	
	get2Dpolymodel(argv[arg], &l, &m, &ndim, &a, &nmodes, &nvar, vardef, &xname);
	if (ndim != 2) {
		error_exit("warpcat: I need a 2D parameter file!\n");
	}

/*	
	for (mode = 0; mode < nmodes; mode++) {
		fprintf(stderr, "%5d %5d %13g %13g\n", l[mode], m[mode], a[0][mode], a[1][mode]);
	}
	exit(0);
*/


        inputcathead = readcathead();                   /* read the cat head */
	getcatipfiletype(&cattype);
 	setcatopfiletype(cattype);
        inputobject = newobject(inputcathead);          /* make the input object */
        connectobjecttocathead(inputobject);            /* obj addresses point back to cathead */
        allocobjectcontents(inputobject);               /* and allocate space for obj data */

        outputcathead = (cathead *) calloc(1, sizeof(cathead));                 /* new cathead */
        copyheaderinfo(outputcathead, inputcathead);                            /* copy header stuff */
        addargscomment(argc, argv, outputcathead);              /* add history */
        copycontentinfo(outputcathead, inputcathead);           /* copy over pre-exisiting object items */
        ritem = newitem("r", NUM_TYPE, 1, 2);                   /* add r object item */
        addobjectitem(ritem, outputcathead);
        writecathead(outputcathead);                            /* and write cathead out */

        outputobject = newobject(outputcathead);                /* make the output object */
        rindex = getobjectitemindex("r", outputobject);		/* get indices for new items */
        inheritcontents(outputobject, inputobject);             /* pre-existing output object addresses point back to inputobject */
        allocitemcontents(ritem, &((outputobject->addrlist)[rindex]), 0);     /* allocate space for new data */
        r = (double *) ((outputobject->addrlist)[rindex]);      /* get local handles for the new items. */
        /* now we get the handles to the input object items we will need */
        x = (double *) ((inputobject->addrlist)[getobjectitemindex("x", inputobject)]);


        while (readobject(inputobject)) {                                       /* big loop */
		rmodel(x, r);
		if (dophitrans) {
			for (i = 0; i < 2; i++) {
				rp[i] = phi[i][0] * r[0] + phi[i][1] * r[1];
			}
			for (i = 0; i < 2; i++) {
				r[i] = rp[i];
			}
		}
		if (dodisplacement) {
			for (i = 0; i < 2; i++) {
				r[i] += d[i];
			}
		}
                writeobject(outputobject);
	}
	exit(0);
}


void	rmodel(double *x, double *r)
{
	int	i, mode;
	double	fmode;

	for (i = 0; i < 2; i++) {
		r[i] = x[i];
	}
	for (mode = 0; mode < nmodes; mode++) {
		fmode = f(l[mode], m[mode], x);
		for (i = 0; i < 2; i++) {
			r[i] += a[i][mode] * fmode;
		}
	}
}
