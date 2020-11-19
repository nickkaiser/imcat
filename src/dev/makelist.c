/*
 * makelist.c
 */

#include <stdio.h>
#include <math.h>

#define usage "usage: makelist nobj [xmin xmax ymin ymax [seed]]\n\
	lay down random objects within a rectangle.\n\n"

double	drand48(), x, y;

main(int argc, char *argv[])
{
	int	seed, nobj;
	double	xmin, xmax, ymin, ymax;

	seed = 1;
	xmin = 0.0;
	xmax = 1.0;
	ymin = 0.0;
	ymax = 1.0;

	if ((argc == 2) && (argv[1][0] = '-')) {
		fprintf(stderr, "%s", usage);
		exit(-1);
	}
	
	switch (argc) {
		case 7:
			sscanf(argv[6], "%d", &seed);
		case 6:
			sscanf(argv[2], "%lf", &xmin);
			sscanf(argv[3], "%lf", &xmax);
			sscanf(argv[4], "%lf", &ymin);
			sscanf(argv[5], "%lf", &ymax);
		case 2:
			if (1 == sscanf(argv[1], "%d", &nobj))
				break;
		case 1:
		default:
			fprintf(stderr, "%s", usage);
			exit(-1);
			break;
	}

	srand48(seed);
	fprintf(stdout, "#        x          y\n");
	while (nobj) {
		fprintf(stdout, "%10.3lf %10.3lf\n", 
			xmin + (xmax - xmin) * drand48(), 
			ymin + (ymax - ymin) * drand48());
		nobj--;
	}
}

