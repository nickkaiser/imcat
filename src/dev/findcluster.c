/*
 * findcluster.c
 */
#include <stdio.h>
#include <math.h>

#define usage "usage: findcluster ncol sigma_1...sigma_ncol [options...]\n\
	searches for a cluster of particles in multi-dimensional space.\n\
		-b dbox		# half-box side d_i = dbox * sigma_i (2.0)\n\
		-w 		# weight the points by f\n\
		-p 		# print the points in the box\n\
	You supply ncol, the number of columns and a smoothing scale\n\
	sigma_i for each dimension.\n\
	For each point we calculate f = sum over neighbours of\n\
	exp(-0.5 sum_i dx_i^2 / sigma_i^2) and find the hottest particle.\n\
	Finally we average coordinate summing over particles\n\
	in the box of half-side d_i = dbox * sigma_i.\n\
	Use -w option to weight points by f.\n\
	Use -p option to print the points in the box (with f-values\n\
	in last column).\n\n"

typedef struct particle {
	double			*x;
	double			f;
	struct	particle	*next;
} particle;

#define	MAX_DIMS 6
#define	D2_MAX	20

main(int argc, char *argv[])
{
	char		line[1024];
	int		arg, printmode, weightbyf, ncol, c, npairs = 0, inbox;
	double		*sigma, dbox, d2, d, fmax, weightsum, *xbar;
	particle	*p, *pbase, *p1, *p2, *pmax;

	/* defaults */
	printmode = 0;
	weightbyf = 0;
	dbox = 2;

	/* parse args */
	if (argc < 3) {
		fprintf(stderr, "%s", usage);
		exit(-1);
	}
	sscanf(argv[1], "%d", &ncol);
	if (ncol > 6 || argc < 2 + ncol) {
		fprintf(stderr, "%s", usage);
		exit(-1);
	}
	sigma = (double *) calloc(ncol, sizeof(double));
	xbar = (double *) calloc(ncol, sizeof(double));
	for (c = 0; c < ncol; c++) {
		sscanf(argv[2 + c], "%lf", &(sigma[c]));
	}
	arg = 2 + c;
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, "%s", usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			case 'b':
				sscanf(argv[arg++], "%lf", &dbox);
				break;
			case 'w':
				weightbyf = 1;
				break;
			case 'p':
				printmode = 1;
				break;
			default:
				fprintf(stderr, "%s", usage);
				exit(-1);
		}
	}

	/* create a NTLL of particles */
	p = NULL;
	while (fgets(line, 1024, stdin)) {
		if (line[0] == '#')
			continue;
		pbase = (particle *) calloc(1, sizeof(particle));
		pbase->x = (double *) calloc(ncol, sizeof(double));
		pbase->next = p;
		p = pbase;
		switch (ncol) {
			case 1:
				sscanf(line, "%lf", p->x);
				break;
			case 2:
				sscanf(line, "%lf %lf", p->x, p->x+1);
				break;
			case 3:
				sscanf(line, "%lf %lf %lf", p->x, p->x+1, p->x+2);
				break;
			case 4:
				sscanf(line, "%lf %lf %lf %lf", p->x, p->x+1, p->x+2, p->x+3);
				break;
			case 5:
				sscanf(line, "%lf %lf %lf %lf %lf", p->x, p->x+1, p->x+2, p->x+3, p->x+4);
				break;
			case 6:
				sscanf(line, "%lf %lf %lf %lf %lf %lf", p->x, p->x+1, p->x+2, p->x+3, p->x+4, p->x+5);
				break;
			default:
				fprintf(stderr, "findcluster: bad ncol\n");
				exit(-1);
				break;
		}
	}
	
	/* now sum over neighbours */
	p1 = pbase;
	while (p1) {
		p2 = pbase;
		while (p2) {
			if (p1 != p2) {
				d2 = 0.0;
				for (c = 0; c < ncol; c++) {
					d = (p2->x[c] - p1->x[c]) / sigma[c];
					d2 += d * d;
				}
				if (d2 < D2_MAX)
					p1->f += exp(-0.5 * d2);
				npairs++;
			}
			p2 = p2->next;
		}
		p1 = p1->next;
	}

	/* find the hottest particle */
	p = pbase;
	fmax = 0.0;
	while (p) {
		if (p->f > fmax) {
			pmax = p;
			fmax = p->f;
		}
		p = p->next;
	}

	/* print the cluster if we're in printmode */
	if (printmode) {
		fprintf(stdout, "# fmax = %lf at\n# ", fmax);
		for (c = 0; c < ncol; c++)
			fprintf(stdout, "%10.3lf ", pmax->x[c]);
		fprintf(stdout, "\n");
	}
	p = pbase;
	weightsum = 0.0;
	while (p) {
		inbox = 1;
		for (c = 0; c < ncol; c++) {
			if (fabs(p->x[c] - pmax->x[c]) > dbox * sigma[c]) {
				inbox = 0;
				break;
			}
		}
		if (inbox) {
			if (printmode) {
				for (c = 0; c < ncol; c++)
					fprintf(stdout, "%10.3lf ", p->x[c]);
				fprintf(stdout, "%10.3lf\n", p->f);
			} else {
				weightsum += (weightbyf ? p->f : 1.0);
				for (c = 0; c < ncol; c++)
					xbar[c] += (weightbyf ? p->f * p->x[c] : p->x[c]);
			}
		}
		p = p->next;
	}
	if (!printmode) {
		for (c = 0; c < ncol; c++) {
			fprintf(stdout, "%13.6lf ", xbar[c] / weightsum);
		}
		fprintf(stdout, "\n");
	}
} 

