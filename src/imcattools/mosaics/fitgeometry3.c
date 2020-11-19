/*
 * fitgeometry3.c
 *
 */


#define usage "\n\
NAME\n\
	fitgeometry3 --- fit layout of a set of images\n\
\n\
SYNOPSIS\n\
	fitgeometry3 [options...]\n\
		-l lmax		# maximum order for polynomial distortion model (1)\n\
		-o outputdir	# directory for the output files (must exist)\n\
		-v		# 'verbose' mode\n\
\n\
DESCRIPTION\n\
	fitgeometry3 is very similar to fitgeometry2.\n\
	The difference is that fitgemetry3 reads a catalogue\n\
	(usually a concatenation of catalogues from a number of images)\n\
	containing, in addition to the pixel position 'x', image number 'i'\n\
	and position measurement variance 's', a particle identifier 'p'\n\
	which is an integer.\n\
\n\
	Output files are put in 'geofit3files/' by default, but you can change\n\
	this with the -o option.\n\
\n\
	With -v option, fitgeometry3 will report how many objects read, size of\n\
	system of equations etc to stderr.\n\
\n\
	See <a href=\"fitgeom.html\">notes on various fitgeometry tools</a>\n\
\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils/modefunc.h"
#include "utils/lu.h"
#include "utils/ipbuff.h"

/* input buffer details */
#define IPRECUNIT	5
#define IPxPOS		0
#define IPiPOS		2
#define IPpPOS		3
#define IPsPOS		4

#define MAXIMAGES	1000
#define MAXPARTICLES	10000

	

main(int argc, char *argv[])
{
	int	arg = 1;		/* arg counter */
	int	nm,			/* number of distortion parameters per image */
		nmodes,			/* total number of distortion parameters */
		np,			/* number of distinct particles */
		nobj,			/* number of objects read from input cat */
		n,			/* total number of equations (n = nmodes + np) */
		lpmax,			/* max order of polynomials */
		*lpa,			/* array of lp-values */
		*mpa,			/* array of mp-values */
		*indx,			/* used by lu routines */
		I,J,i,j,l,m,lp,mp,p,q,r,/* misc indices */
		iprecsize,		/* input record size*/
		nim,			/* number of images */
		*ndet, *refdet,		/* number of detections; was it detected in reference cat */
		*imind,	*imlab,		/* index array and label array for images */
		verbose,		/* output useful stats if(verbose) */
		*pind, *plab;		/* particle number array and 'inverse' or label array */
	double	*a[2],			/* transformation parameter array */		
		**A, *B[2], det,	/* linear algebra matrix A[I][J], rhs vector B[2][J], determinant*/
		*x,			/* measured coordinates */
		s,			/* measurement variance */
		*r0, *r1,		/* components of the solution */
		**buffp;		/* input buffer */
	FILE	*lcpipe; 		/* input catalogue */
	char	opfilename[1024],	/* parameter file */
		*opfilebase,		/* prefix */
		*vardef[1],		/* write2Dpolymodel needs array of vardefs as argument */
		pipestring[512],	/* temporary string */
		tmpstring[512];		/* temporary string */

	/* defaults */
	lpmax = 1;
	vardef[0] = "1 2 a";
	opfilebase = "geofit3files";
	verbose = 0;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
			case 'l':
				if (argc - arg < 1) {
					error_exit(usage);
				}
				if ((1 != sscanf(argv[arg++], "%d", &lpmax))) {
						error_exit(usage);
				}
				break;
			case 'o':
				if (argc - arg < 1) {
					error_exit(usage);
				}
				opfilebase = argv[arg++];
				break;
			case 'v':
				verbose = 1;
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}

	/* check to see that we can write to opfilebase */
	sprintf(tmpstring, "test -d %s", opfilebase);
	if (system(tmpstring)) {
		fprintf(stderr, "fitgeometry3: you need to create the output file directory \"%s/\" first.\n", opfilebase);
		exit(-1);
	}

	/* copy args to modefunc argstring */
	modefunc_addargcomment(argc, argv);

	/* compute number of coefficients per chip*/
	nm = 0;
	for (lp = 0; lp <= lpmax; lp++) {
		nm += (lp + 1);
	}
	/* and set up arrays of lp, mp values */
	lpa = (int *) calloc(nm, sizeof(int));
	mpa = (int *) calloc(nm, sizeof(int));
	i = 0;
	for (lp = 0; lp <= lpmax; lp++) {
		for (mp = 0; mp <= lp; mp++) {
			lpa[i] = lp;
			mpa[i] = mp;
			i++;
		}			
	}

	/* open lc pipe for input */
	iprecsize = IPRECUNIT;
        if (!(lcpipe = popen("lc -b -o x i p s", "r"))) {
                fprintf(stderr, "fitgeometry3: unable to open pipe 'lc -b -o x i p s' for input\n");
                exit(-1);
        }
	/* and read into an expandable buffer */
	buffp = readdoublebuff(IPRECUNIT, lcpipe, &nobj);
	pclose(lcpipe);

	/* make the array of image indices imind[] and compute number of images*/
	imind = (int *) calloc(MAXIMAGES, sizeof(int));
	imlab = (int *) calloc(MAXIMAGES, sizeof(int));
	for (p = 0; p < nobj; p++) {
		i = (int) buffp[p][IPiPOS];
		if (i >= MAXIMAGES) {
			error_exit("fitgeometry3.c: image number too large\n");
		}
		imind[i] = 1;
	}
	imind[0] = imlab[0] = 0;
	i = 1;
	for (j = 1; j < MAXIMAGES; j++) {
		if (imind[j]) {
			imind[j] = i;
			imlab[i] = j;
			i++;
		}
	}	
	nim = i;	/* total number of images (including reference image) */

	/* make the array of particle numbers and the inverse array */
	pind = (int *) calloc(MAXPARTICLES, sizeof(int));
	plab = (int *) calloc(MAXPARTICLES, sizeof(int));
	for (p = 0; p < nobj; p++) {
		q = (int) buffp[p][IPpPOS];
		if (q >= MAXPARTICLES) {
			error_exit("fitgeometry3.c: particle number too large\n");
		}
		pind[q] = 1;
	}
	p = 0;
	for (q = 0; q < MAXPARTICLES; q++) {
		if (pind[q]) {
			pind[q] = p;
			plab[p] = q;
			p++;
		}
	}
	np = p;		/* total number of distinct particles */

	/* now we can allocate and compute ndet, refdet */
	ndet = (int *) calloc(np, sizeof(int));
	refdet = (int *) calloc(np, sizeof(int));
	for (q = 0; q < nobj; q++) {
		p = pind[(int) buffp[q][IPpPOS]];
		ndet[p]++;
		i = imind[(int) buffp[q][IPiPOS]];
		if (!i) {
			refdet[p] = 1;
		}
	}
	
	/* now we can allocate arrays for linear algebra arrays */
	nmodes = nm * (nim - 1);
	n = nmodes + np;
	
	/* output useful stats */
	if (verbose) {
		fprintf(stderr, "fitgeometry3: %d objects read\n", nobj);
		fprintf(stderr, "fitgeometry3: %d distinct particles\n", np);
		fprintf(stderr, "fitgeometry3: %d distinct images\n", nim);
		fprintf(stderr, "fitgeometry3: %d linear equations\n", n);		
	}

       	A = (double **) calloc(n, sizeof(double *));
 	for (I = 0; I < n; I++) {
                A[I] = (double *) calloc(n, sizeof(double));
       	}	
	for (q = 0; q < 2; q++) {
        	B[q] = (double *) calloc(n, sizeof(double));
 	}

	/* accumulate  A-matrix, B-vector */
	for (q = 0; q < nobj; q++) {
		x = buffp[q] + IPxPOS;
		i = imind[(int) buffp[q][IPiPOS]];
		p = pind[(int) buffp[q][IPpPOS]];
		s = buffp[q][IPsPOS];
		/* A-matrix */
		if (i) {
		/* upper left */
			for (l = 0; l < nm; l++) {
				for (m = 0; m < nm; m++) {
					I = (i - 1) * nm + l;
					J = (i - 1) * nm + m;
					A[J][I] += f(lpa[l], mpa[l], x) * f(lpa[m], mpa[m], x) / s;
				}
			}
			/* upper right */
			for (m = 0; m < nm; m++) {
				J = (i - 1) * nm + m;
				I = nmodes + p;
				if (A[J][I] != 0.0) {
					error_exit("fitgeometry3: object appears multiply on single image!\n");
				}
				A[J][I] = - f(lpa[m], mpa[m], x) / s;
			}
			/* lower left */
			for (l = 0; l < nm; l++) {
				I = (i - 1) * nm + l;
				J = nmodes + p;
				if (A[J][I] != 0.0) {
					error_exit("fitgeometry3: object appears multiply on single image!\n");
				}
				A[J][I] = f(lpa[l], mpa[l], x) / s;
			}
		}
		/* bottom right */
		I = J = nmodes + p;
		A[J][I] -= 1.0 / s;
		/* B-vectors */
		for (r = 0; r < 2; r++) {
			/* upper portion */
			if (i) {
				for (m = 0; m < nm; m++) {
					J = (i - 1) * nm + m;
					B[r][J] -= f(lpa[m], mpa[m], x) * x[r] / s;
				}
			}
			/* lower portion */
			J = nmodes + p;
			B[r][J] -= x[r] / s;
		}
	}



	
  	/* solve the linear equations */
	indx = (int *) calloc(n, sizeof(int));
	myludcmp(A, n, indx, &det);
	for (q = 0; q < 2; q++) {
		mylubksb(A, n, indx, B[q]);
	}

	/* extract the model coefficients and write the  parameter files */
	for (i = 1; i < nim; i++) {
		a[0] = B[0] + (i - 1) * nm;
		a[1] = B[1] + (i - 1) * nm;
		sprintf(opfilename, "%s/%d%s", opfilebase, imlab[i], ".par");
		write2Dpolymodel(opfilename, nm, lpa, mpa, 2, a, 1, vardef, "x");
	}

	/* generate the grand reference catalogue */
	r0 = B[0] + nmodes;
	r1 = B[1] + nmodes;
	sprintf(tmpstring, "lc -C -b -x -a 'generated by fitgeometry3' -N '1 2 r' -n p -n ndet -n refdet > %s/ref.cat", opfilebase);
	lcpipe = popen(tmpstring, "w");
	if (!lcpipe) {
		error_exit("fitgeometry3: unable to open lcpipe for output\n");
	}
	for (p = 0; p < np; p++) {
		fprintf(lcpipe, "%13.8lg %13.8lg %d %d %d\n", r0[p], r1[p], plab[p], ndet[p], refdet[p]);
	}	
	pclose(lcpipe);
	/* and generate the output catalogues for each image */
	for (i = 0; i < nim; i++) {
		sprintf(opfilename, "%s/%d%s", opfilebase, imlab[i], ".par");
		if (i == 0) {
			 sprintf(pipestring, "| lc -b +all 'r = %%x' > ");
		} else {
			sprintf(pipestring, "| warpcat %s > ", opfilename);
		}
		sprintf(opfilename, "%s/%d%s", opfilebase, imlab[i], ".cat");
		sprintf(tmpstring, "%s %s %s",
			"lc -C -b -x -a 'generated by fitgeometry3' -N '1 2 x' -n i -n p -n s -N '1 2 rref' -n ndet -n refdet ",
			pipestring, opfilename);
		lcpipe = popen(tmpstring, "w");
		if (!lcpipe) {
			error_exit("fitgeometry2: unable to open lcpipe for output\n");
		}
		for (q = 0; q < nobj; q++) {
			j = imind[(int) buffp[q][IPiPOS]];
			if (j == i) {
				x = buffp[q] + IPxPOS;
				p = pind[(int) buffp[q][IPpPOS]];
				s = buffp[q][IPsPOS];
				fprintf(lcpipe, "%13.8lg %13.8lg %1d %1d %13.8lg %13.8lg %13.8lg %d %d\n",
						x[0], x[1], imlab[i], plab[p], s, r0[p], r1[p], ndet[p], refdet[p]);
			}
		}
		pclose(lcpipe);
	}
	

	/* all done */
	exit(0);
}
