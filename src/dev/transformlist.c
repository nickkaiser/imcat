/*
 * transformlist.c
 */
#include <stdio.h>
#include <math.h>

#define usage "usage: transformlist a phi dx dy\n\
	applies transformation\n\
		x0 = a (x cos phi - y sin phi) + dx\n\
		y0 = a (x sin phi + y cos phi) + dy\n\n"

main(int argc, char *argv[])
{
	double	x, y, x0, y0, a, phi, dx, dy;
	char	line[1024];

	if (argc != 5) {
		fprintf(stderr, "%s", usage);
		exit(-1);
	}
	sscanf(argv[1], "%lf", &a);
	sscanf(argv[2], "%lf", &phi);
	sscanf(argv[3], "%lf", &dx);
	sscanf(argv[4], "%lf", &dy);

	while (fgets(line, 1024, stdin)) {
		if (line[0] == '#')
			continue;
		sscanf(line, "%lf %lf ", &x, &y);
		x0 = a * (x * cos(phi) - y * sin(phi)) + dx;
		y0 = a * (x * sin(phi) + y * cos(phi)) + dy;
		fprintf(stdout, "%10.3lf %10.3lf\n", x0, y0);
	}
}

