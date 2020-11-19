/*
 * fit2cats.c
 */


#define usage "\n\
NAME\n\
	fit2cats --- fit for transformation/distortion for a pair of cats\n\
\n\
SYNOPSIS\n\
	fit2cats \n\
		-l lmax		# max order for distortion polynomials (1)\n\
		-o x0 y0	# origin for spatial coordinates (0,0)\n\
\n\
DESCRIPTION\n\
	'fit2cats' reads from stdin the result from 'merge2cats' of\n\
	merging a pairs of cats, and which must contain entries for\n\
	a pair of spatial coordinate vectors 'x[2][2]'.\n\
	It then fits a model in which r = x[1] is related to (x,y) = x[0]\n\
		r = x[] + sum_{l=0}^{l_max} sum_{m=0}^l a_{lm} x^{l-m} y^m\n\
	With the default l_max = 1, this is just a straight linear\n\
	transformation, but by going to higher order one can include\n\
	distortions of the telescope optics.\n\
	With -o option, these mode functios become polynomials in\n\
	position relative to specified spatial origin.\n\
	Transformation parameters a_{lm} are written to stdout in form\n\
	to be read by warpimage or warpcat.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils/modefunc.h"

#define SCALE 1.0

/* nprs-length vectors of data */
double 		***x;


int     readmergedcat(void);
	
main(int argc, char *argv[])
{
	/* number of pairs, exposures, etc */
	static	int	nprs, nmodes, lmax, *ll, *mm, *indx;
	/* model parameters */
	static	double	*a[2];
	/* matrix stuff */
	static	double	**A, *B, det;
	/* arg counter */
	int	arg = 1;
	int	i, j, l, m, n, M, ipr;
	double	x0[2];
	char	*vardef[1];

	/* defaults */
	lmax = 1;
	x0[0] = x0[1] = 0.0;
	vardef[0] = "1 2 a";

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
				if ((1 != sscanf(argv[arg++], "%d", &lmax))) {
						error_exit(usage);
				}
				break;
			case 'o':
				if (argc - arg < 2) {
					error_exit(usage);
				}
				if (1 != sscanf(argv[arg++], "%lf", &(x0[0])) ||
					1 != sscanf(argv[arg++], "%lf", &(x0[1]))) {
						error_exit(usage);
				}
				setorigin(x0);
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}	

	/* copy args to modefunc argstring */
	modefunc_addargcomment(argc, argv);

	/* read the merged catalogue */
	nprs = readmergedcat();
	fprintf(stderr, "# read %d pairs\n", nprs);

	/* compute number of coefficients */
	nmodes = 0;
	for (l = 0; l <= lmax; l++) {
		nmodes += (l + 1);
	}

	/* allocate arrays for linear algebra */
        B = (double *) calloc(nmodes, sizeof(double));
        A = (double **) calloc(nmodes, sizeof(double *));
        for (m = 0; m < nmodes; m++) {
                A[m] = (double *) calloc(nmodes, sizeof(double));
        }
        indx = (int *) calloc(nmodes, sizeof(int));

	/* allocate space for the model parameters */
	for (i = 0; i < 2; i++) {
		a[i] = (double *) calloc(nmodes, sizeof(double));
	}

	/* set up arrays of l, m values */
	ll = (int *) calloc(nmodes, sizeof(int));
	mm = (int *) calloc(nmodes, sizeof(int));
	M = 0;
	for (l = 0; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			ll[M] = l;
			mm[M++] = m;
		}			
	}	

	/* solve for a[i], i=0,1 */
	for (i = 0; i < 2; i++) {
		/* zero arrays */
		for (m = 0; m < nmodes; m++) {
			B[m] = 0.0;
			for (n = 0; n < nmodes; n++) {
				A[m][n] = 0.0;
			}
		}
		/* accumulate arrays */
        	for (ipr = 0; ipr < nprs; ipr++) {
			for (m = 0; m < nmodes; m++) {
				B[m] += (x[ipr][1][i] - x[ipr][0][i]) * f(ll[m], mm[m], x[ipr][0]);
				for (n = 0; n < nmodes; n++) {
					A[m][n] +=  f(ll[m], mm[m], x[ipr][0]) * f(ll[n], mm[n], x[ipr][0]);
				}
			}
		}
		/* solve the linear equations */
		myludcmp(A, nmodes, indx, &det);
		mylubksb(A, nmodes, indx, B);
		/* extract the model coefficients */
		for (m = 0; m < nmodes; m++) {
			a[i][m] = B[m];
			/* subtract origin shift 
			if (ll[m] == 0) {
				a[i][m] -= x0[i];
			}*/
			/* subtract identity
			if ((ll[m] == 1) && (mm[m] == i)) {
				a[i][m] -= 1.0;
			} */
		}
	}

	write2Dpolymodel(NULL, nmodes, ll, mm, 2, a, 1, vardef, "x");
	exit(0);
}





/* number of objects per input buffer */
#define BUFF_SIZE 1000
/* number of fields */
#define N_FIELDS 4
/* max number of input buffers */
#define MAX_BUFFS 10000
int     readmergedcat(void)
{
        FILE    *lcpipe;
        int     i, j, ibuff, nbuff, nprs, ipr, maxprs = BUFF_SIZE, nprsinbuff;
        double  *buff[MAX_BUFFS];

        /* open pipe */
        if (!(lcpipe = popen("lc -b -o x", "r"))) {
                fprintf(stderr, "fit2cats: readmergedcat: unable to open lc-pipe for input\n");
                exit(-1);
        }

        /* read catalogue into array of buffers */
        nbuff = 0;
        nprs = nprsinbuff = 0;
        if (!(buff[nbuff] = (double *) calloc(N_FIELDS * BUFF_SIZE, sizeof(double)))) {
                error_exit("fit2cats: readmergedcat: failed to allocate input buffer\n");
        }
        while (fread(buff[nbuff] + nprsinbuff * N_FIELDS, sizeof(double), 4, lcpipe)) {
                nprsinbuff++;
                nprs++;
                if (nprsinbuff == BUFF_SIZE) {  /* need a new buffer */
                        if (nbuff++ >= MAX_BUFFS) {
                                error_exit("fit2cats: readmergedcat: too many objects in input cat\n");
                        }
                        if (!(buff[nbuff] = (double *) calloc(N_FIELDS * BUFF_SIZE, sizeof(double)))) {
                                error_exit("fit2cats: readmergedcat: failed to allocate input buffer\n");
                        }
                        nprsinbuff = 0;
                }
        }
        nbuff++;

        /* create data arrays */
        x = (double ***) calloc(nprs, sizeof(double **));
        for (ipr = 0; ipr < nprs; ipr++) {
                x[ipr] = (double **) calloc(2, sizeof(double *));
                 for (j = 0; j < 2; j++) {
                        x[ipr][j] = (double *) calloc(2, sizeof(double));
                }
        }

        /* copy data from input buffers */
        nprs = 0;
        for (ibuff = 0; ibuff < nbuff; ibuff++) {
                if (ibuff == (nbuff - 1)) {
                        maxprs = nprsinbuff;
                } else {
                        maxprs = BUFF_SIZE;
                }
                for (ipr = 0; ipr < maxprs; ipr++) {
                        x[nprs][0][0] = buff[ibuff][ipr * N_FIELDS];
                        x[nprs][0][1] = buff[ibuff][ipr * N_FIELDS + 1];
                        x[nprs][1][0] = buff[ibuff][ipr * N_FIELDS + 2];
                        x[nprs][1][1] = buff[ibuff][ipr * N_FIELDS + 3];
                        nprs++;
                }
                free(buff[ibuff]);
        }
        return (nprs);
}  

