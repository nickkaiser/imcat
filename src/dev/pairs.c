/*
 * pairs.c
 *
 */


#include <stdio.h>
#include <math.h>

#define usage "usage: pairs\n\
	reads a list of x,y,m(agnitude) and outputs list of pairs of points.\n\
	Outputs x1,y1,m1,dr,phi,dm\n\
	Where x1,y1,m1 are pos,mag of brightest object\n\
	and dr, phi, dm are |r2-r1|, atan2(y2-y1,x2-x1), m2-m1.\n\n"

#define		WHOPPER		1.e10
#define 	MAX_OBJECTS	10000

#define	dist(i,j)	sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]))

double	x[MAX_OBJECTS], y[MAX_OBJECTS], m[MAX_OBJECTS];

main(int argc, char *argv[])
{
	char	line[1024];
	int	n, i1, i2;
	double	dr, phi;

	if (argc != 1) {
		fprintf(stderr, "%s", usage);
		exit(-1);
	}

	/* read the objects */
	n = 0;
	while (fgets(line, 1024, stdin)) {
		if (line[0] == '#')
			continue;
		sscanf(line, "%lf %lf %lf", x + n, y + n, m + n);
		n++;
	}

	/* find the pairs */
	fprintf(stdout, "#       x1         y1         m1         dr        phi         dm\n");
	for (i1 = 0; i1 < n; i1++) {
		for (i2 = 0; i2 < n; i2++) {
			if (i1 == i2 || m[i1] > m[i2])
				continue;
			dr = dist(i1,i2);
			phi=atan2(y[i2]-y[i1],x[i2]-x[i1]);
			fprintf(stdout, "%10.3lf %10.3lf %10.3lf %10.3lf %10.5lf %10.3lf\n",
				x[i1], y[i1], m[i1], dr, phi, m[i2] - m[i1]);
		}
	}
}
