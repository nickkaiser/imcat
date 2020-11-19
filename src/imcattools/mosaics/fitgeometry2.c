

/*
 * fitgeometry2.c
 *
 *  In this variant we read in the result of merging a set of 'planes' of images,
 *  the first of which should be a 'reference' image (with plane number = 0).
 *  We then solve for a set of transformations for all the other images and
 *  for the positions of the N distinct objects which appear in the input
 *  merged catalogue.  
 *
 *  Let there be nm modes for for each image: i.e. coordinates in the
 *  reference system for the p'th object in terms of its measured
 *  position on the i'th image (if present) is given by
 *     		r_p = x_pi + sum_m a_im f_m(x_pi) + e_pi
 *  whereas for the reference image
 *     		r_p = x_p0 + e_p0
 *  where the e_pi is the error in the p'th object's position on the ith frame
 *  with <e_pi^2> = s_pi.
 *
 *  The chi-squared we need to minimise is
 * 
 *	chi^2 = sum_p [ sum_i E_pi (x_pi - sum_l a_il f_l(x_pi) - r_p)^2 / s_pi^2
 *			+ E_p0 (x_p0 - r_p)^2 / s_po]
 *  where E_pi = 1 if object p is found in catalgoue i.
 *
 *  Minimisation gives
 *	nm * (nimages - 1) + np
 *  linear equations for the distortion mode amplitudes a_il and
 *  the np object locations.  The result is a set of transformation parameter
 *  files and locations for a 'reference catalogue' which can be then used
 *  to warp images.
 *
 *  This gives the following set of linear equations for the a_il, r_p's:
 *
 *	A x = B = 
 *	= -   - 
 *
 *  -----------------------------------  --------     --------
 * |            :                      ||        |   |        |
 * | A_(jm)(il) :       A_(jm)(p)      || a_(il) |   | B_(jm) |
 * |            :                      ||        |   |        |
 * |...................................||........| = |........|
 * |            :                      ||        |   |        |
 * |            :                      ||        |   |        |
 * |            :                      ||        |   |        |
 * |            :                      ||        |   |        |
 * | A_(p')(il) :       A_(p')(p)      ||  r_(p) |   | B_(p') |
 * |            :                      ||        |   |        |
 * |            :                      ||        |   |        |
 * |            :                      ||        |   |        |
 * |            :                      ||        |   |        |
 *  -----------------------------------  --------     --------
 *
 * where
 *
 *      A_(jm)(il) = delta_ij sum_p E_pi f_l(x_pi) f_m(x_i) / s_pi
 *
 *      A_(jm)(p)  = - E_pj f_m(x_pj) / s_pj
 *
 *      A_(p')(il) = E_p'i f_l(x_p'i) / s_p'i
 *
 *      A_(p')(p)  = - delta_(p)(p') sum_i=0 1 / s_pi	
 *
 *      B_(jm)     = - sum_p f_m(x_pj) x_pj / s_pj
 *
 *      B_(p')     = - sum_i=0 E_p'i x_p'i / s_p'i
 *
 */


#define usage "\n\
NAME\n\
	fitgeometry2 --- fit layout of a set of images\n\
\n\
SYNOPSIS\n\
	fitgeometry2 np [options...]\n\
		-l lmax		# maximum order for polynomial distortion model (1)\n\
		-o outputdir	# directory for the output files (must exist)\n\
\n\
DESCRIPTION\n\
	fitgeometry2 reads from stdin the result of merging (using 'mergecats')\n\
	a set of np 'planes' of catalogues and solves for the location\n\
	in some 'reference frame coordinates' of objects on the catalogue\n\
	and also a set of parameters describing the distorted mapping between\n\
	these reference frame coordinates and pixel coordinates on the\n\
	images from which the catalogues were derived.\n\
\n\
	Fitgeometry2 was written to solve the following problem: We have a set\n\
	of 'data' images from a mosiac camera (which have some uncertain\n\
	layout of the chips on the detector frame - which need not be static - and\n\
	also suffer from telescope field distortion and possibly atmospheric refraction).\n\
	From each these images one can extract typically ~100 stars whose positions\n\
	have a precision of a small fraction <~ 1/10 of a pixel. These allow one to\n\
	determine a set of polynomial mappings (with one image taken to define\n\
	the reference coordinate system) which map these images onto one\n\
	another to a very high precision.  However, this procedure does not\n\
	remove telescope or atmospheric field distortion and tends to be\n\
	unstable to introducing further artificial field distortion.  To avoid this\n\
	we incorporate in the fitting a catalogue which derives from the digital\n\
	sky survey image for example, and for which a good 'plate solution'\n\
	already exists.  Providing this catalogue as the 'reference catalogue'\n\
	solves our problem and provides one with a mapping from data pixel\n\
	coordinates to the `world coordinate system'.  A complication of this\n\
	procedure is that the reference catalogue positions tend to be relatively\n\
	imprecise, and it is necessary to incorporate this information in the\n\
	fitting (as weight factors), so the input catalogues most contain both\n\
	the measured position and an estimate of the precision.\n\
\n\
	The catalogues to be merged must contain at least the following items:\n\
\n\
		x[2]	# spatial coordinate\n\
		i	# unique image ID (i = 0 for reference image)\n\
		I	# I = 1 -- used internally\n\
		s	# sigma^2 = position measurement variance   \n\
\n\
	but will usually contain additional information such as an approximate\n\
	'sky coordinate' used by mergecats to link the objects.\n\
	The result of such a merging is a set of unique particles and their\n\
	measured positions in each and all of the images in which they appear.\n\
\n\
	The 'planes' which get merged could be just the set of all images.  However,\n\
	in the case of a mosaic camera an object can only be detected on at most\n\
	one chip per exposure, so the merged catalogue will be very sparse.  It is\n\
	more efficient to group sets of such 'known to be mutually non-overlapping'\n\
	catalogues into planes before merging.  Since one cannot then infer\n\
	the 'image number' of a catalogue from its 'plane number' we include the image\n\
	identifier 'i' in the catalogue.\n\
\n\
	Fitgeometry2 assumes that coordinates x_pi of the p'th of np stars on the i'th\n\
	image is related to the coordinates r_p in the reference frame by:\n\
\n\
     		r_p = x_pi + sum_m a_im f_m(x_pi) + e_pi\n\
\n\
	whereas for the reference image (i = 0)\n\
\n\
     		r_p = x_p0 + e_p0\n\
\n\
	Here the f_m(x) are a set of nm polynomial mode functions\n\
	and a_im denotes the amplitudes of these functions in the distortion\n\
	of the i'th image.  The e_pi represent the uncertainty in the position.\n\
	It finds the set of parameters a_il and positions r_p which minimise\n\
	the 'chi-squared' function:\n\
\n\
	chi^2 = sum_p [ sum_i (x_pi - sum_l a_il f_l(x_pi) - r_p)^2 / s_pi^2\n\
			+ (x_p0 - r_p)^2 / s_p0]\n\
\n\
	Since this is quadratic in the a_il, r_p the result is a set of\n\
	nm * (nimages - 1) + np linear equations for each of the 2 spatial\n\
	coordinate components.\n\
	The summation over 'i' here does not include the reference\n\
	catalogue (i=0). For small distortions this is equivalent to the\n\
	maximum likelihood solution if we assume that the position errors are\n\
	gaussian distributed.\n\
\n\
	The result is:\n\
	1) A set of parameter files ouputdir/i.par, where i is the image number.\n\
	By default ouputdir = \"geofit2dir\", but you can change this\n\
	with the -o flag.\n\
	2) A set of catalogues ouputdir/i.cat which contains the\n\
	x, i, I, s values for the points detected on the i'th image, as well\n\
	as the following items:\n\
		rref		# the reference catalogue solution (rref = r_p)\n\
		p		# the particle number\n\
		ndet		# the number of images in which particle p was detected\n\
		refdet		# refdet = 1 if particle was detected in the reference cat\n\
		r		# obtained by applying the approximate polynomial model   \n\
	3) A reference catalogue ouputdir/ref.cat containing the solutions r_p.\n\
\n\
        See also <a href=\"fitgeom.html\">notes on various fitgeometry tools</a>\n\
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
#define IPIPOS		3
#define IPsPOS		4

#define MAXIMAGES	1000

main(int argc, char *argv[])
{
	int	arg = 1;		/* arg counter */
	int	nm,			/* number of distortion parameters per image */
		nmodes,			/* total number of distortion parameters */
		np,			/* number of objects */
		n,			/* total number of equations (n = nmodes + np) */
		lpmax,			/* max order of polynomials */
		*lpa,			/* array of lp-values */
		*mpa,			/* array of mp-values */
		*indx,			/* used by lu routines */
		I,J,i,j,l,m,lp,mp,q,	/* misc indices */
		p, pp, e, e1, e2,	/* misc indices */
		iprecsize,		/* input record size*/
		**im,			/* image identifier im[p][e] */
		**exist, ex,		/* particle p exists in image i => exist[p][i] = 1 */
		nplanes, nim,		/* number of exposure planes, images */
		*imind,	*imlab,		/* index array and label array for images */
		ndet;			/* number of times particle was detected */
	double	*a[2],			/* transformation parameter array */		
		**A, *B[2], det,	/* linear algebra matrix A[I][J], rhs vector B[2][J], determinant*/
		***x,			/* measured coordinates x[p][i][2] */
		**s,			/* measurement variance s[p][i] */
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
	opfilebase = "geofit2files";

	/* parse args */
	if (argc < 2) {
		error_exit(usage);
	} else {
		if (argv[arg][0] == '-') {
			error_exit(usage);
		}
	}
	if (1 != sscanf(argv[arg++], "%d", &nplanes)) {
		error_exit(usage);
	}
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
			case 'u':
			default:
				error_exit(usage);
		}
	}

	/* check to see that we can write to opfilebase */
	sprintf(tmpstring, "test -d %s", opfilebase);
	if (system(tmpstring)) {
		fprintf(stderr, "fitgeometry2: you need to create the output file directory \"%s/\" first.\n", opfilebase);
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
	iprecsize = IPRECUNIT * nplanes;
        if (!(lcpipe = popen("lc -b -o x i I s", "r"))) {
                fprintf(stderr, "fitgeometry2: unable to open pipe 'lc -b -o x i I s' for input\n");
                exit(-1);
        }
	/* and read into an expandable buffer */
	buffp = readdoublebuff(nplanes * IPRECUNIT, lcpipe, &np);
	pclose(lcpipe);

	/* make the array of image indices imind[] */
	imind = (int *) calloc(MAXIMAGES, sizeof(int));
	imlab = (int *) calloc(MAXIMAGES, sizeof(int));
	for (p = 0; p < np; p++) {
		for (e = 0; e < nplanes; e++) {
			ex = (int) buffp[p][IPIPOS * nplanes + e];
			if (ex) {
				i = (int) *(buffp[p] + IPiPOS * nplanes + e);
				if (i < MAXIMAGES) {
					imind[i] = 1;
				} else {
					error_exit("fitgeometry2: image number exceeds MAXIMAGES\n");
				}
			}
		}
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

	/* allocate arrays x[p][i][2], exist[p][i], s[p][i] */
	x =     (double ***) calloc(np, sizeof(double **));
	exist = (int **) calloc(np, sizeof(int *));
	s =     (double **) calloc(np, sizeof(double *));
	for (p = 0; p < np; p++) {
		x[p]      = (double **) calloc(nim, sizeof(double *));
		for (i = 0; i < nim; i++) {
			x[p][i] = (double *) calloc(2, sizeof(double));
		}
		exist[p]  = (int *) calloc(nim, sizeof(int));
		s[p]      = (double *) calloc(nim, sizeof(double));
		/* copy data */
		for (e = 0; e < nplanes; e++) {
			ex = (int) *(buffp[p] + IPIPOS * nplanes + e);
			if (ex) {
				i = imind[(int) *(buffp[p] + IPiPOS * nplanes + e)];
				for (q = 0; q < 2; q++) {
					x[p][i][q] = *(buffp[p] + IPxPOS * nplanes + 2 * e + q);
				}
				exist[p][i] = ex;
				s[p][i] = *(buffp[p] + IPsPOS * nplanes + e);
			}
		}		
	}
	
	/* now we can allocate arrays for linear algebra arrays */
	nmodes = nm * (nim - 1);
	n = nmodes + np;
       	A = (double **) calloc(n, sizeof(double *));
 	for (I = 0; I < n; I++) {
                A[I] = (double *) calloc(n, sizeof(double));
       	}	
	for (q = 0; q < 2; q++) {
        	B[q] = (double *) calloc(n, sizeof(double));
 	}

	/* compute A-matrix */
	/* upper left */
	for (i = 1; i < nim; i++) {
		for (l = 0; l < nm; l++) {
			I = (i - 1) * nm + l;
			for (m = 0; m < nm; m++) {
				J = (i - 1) * nm + m;
				for (p = 0; p < np; p++) {
					if (exist[p][i]) {
						A[J][I] += f(lpa[l], mpa[l], x[p][i]) * f(lpa[m], mpa[m], x[p][i]) / s[p][i];
					}
				}
			}
		}
	}
	/* upper right */
	for (j = 1; j < nim; j++) {
		for (m = 0; m < nm; m++) {
			J = (j - 1) * nm + m;
			for (p = 0; p < np; p++) {
				I = nmodes + p;
				if (exist[p][j]) {
					A[J][I] -= f(lpa[m], mpa[m], x[p][j]) / s[p][j];
				}
			}
		}
	}
	/* lower left */
	for (i = 1; i < nim; i++) {
		for (l = 0; l < nm; l++) {
			I = (i - 1) * nm + l;
			for (p = 0; p < np; p++) {
				J = nmodes + p;
				if (exist[p][i]) {
					A[J][I] += f(lpa[l], mpa[l], x[p][i]) / s[p][i];
				}
			}
		}
	}
	/* bottom right */
	for (p = 0; p < np; p++) {
		I = J = nmodes + p;
		for (i = 0; i < nim; i++) {
			if (exist[p][i]) {
				A[J][I] -= 1.0 / s[p][i];
			}
		}
	}

	/* compute B-vectors */
	/* upper portion */
	for (q = 0; q < 2; q++) {
		for (j = 1; j < nim; j++) {
			for (m = 0; m < nm; m++) {
				J = (j - 1) * nm + m;
				for (p = 0; p < np; p++) {
					if (exist[p][j]) {
						B[q][J] -= f(lpa[m], mpa[m], x[p][j]) * x[p][j][q] / s[p][j];
					}
				}
			}
		}
	}
	/* lower portion */
	for (q = 0; q < 2; q++) {
		for (p = 0; p < np; p++) {
			J = nmodes + p;
			for (i = 0; i < nim; i++) {
				if (exist[p][i]) {
					B[q][J] -= x[p][i][q] / s[p][i];
				}
			}
		}
	}
	
  	/* solve the linear equations */
	fprintf(stderr, "# fitgeometry2 : solving %d simulataneous linear equations! - please wait....\n", n);
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
	/* and generate the output catalogues for each image */
	r0 = B[0] + nmodes;
	r1 = B[1] + nmodes;
	for (i = 0; i < nim; i++) {
		sprintf(opfilename, "%s/%d%s", opfilebase, imlab[i], ".par");
		if (i == 0) {
			 sprintf(pipestring, "| lc +all 'r = %%x' > ");
		} else {
			sprintf(pipestring, "| warpcat %s > ", opfilename);
		}
		sprintf(opfilename, "%s/%d%s", opfilebase, imlab[i], ".cat");
		sprintf(tmpstring, "%s %s %s",
			"lc -C -x -a 'generated by fitgeometry2' -N '1 2 x' -n i -n I -n s -N '1 2 rref' -n p -n ndet -n refdet ",
			pipestring, opfilename);
		lcpipe = popen(tmpstring, "w");
		if (!lcpipe) {
			error_exit("fitgeometry2: unable to open lcpipe for output\n");
		}
		for (p = 0; p < np; p++) {
			if (exist[p][i]) {
				ndet = 0;
				for (j = 0; j < nim; j++) {
					ndet += exist[p][j];
				}
				fprintf(lcpipe, "%13.8lg %13.8lg %1d %1d %13.8lg %13.8lg %13.8lg %d %d %d\n",
						x[p][i][0], x[p][i][1], i, exist[p][i], s[p][i], r0[p], r1[p], p, ndet, exist[p][0]);
			}
		}
		pclose(lcpipe);
	}
	/* finally, generate the grand reference catalogue */
	sprintf(tmpstring, "lc -C -x -a 'generated by fitgeometry2' -N '1 2 r' -n p -n ndet -n refdet > %s/ref.cat", opfilebase);
	lcpipe = popen(tmpstring, "w");
	if (!lcpipe) {
		error_exit("fitgeometry2: unable to open lcpipe for output\n");
	}
	for (p = 0; p < np; p++) {
		ndet = 0;
		for (j = 0; j < nim; j++) {
			ndet += exist[p][j];
		}
		fprintf(lcpipe, "%13.8lg %13.8lg %d %d %d\n", r0[p], r1[p], p, ndet, exist[p][0]);
	}	
	pclose(lcpipe);
	

	/* all done */
	exit(0);
}
