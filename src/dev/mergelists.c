/*
 * mergelists.c
 */

#include <stdio.h>
#include <math.h>



#define usage "usage: mergelists list1 list2 [options...]\n\
	reads two lists containing (x1, y2) and (x2, y2) particle coords.\n\
	Construct a N x N checkerboard containing null\n\
	terminated linked lists of particles from list1\n\
	and then use this to find nearest neighbour of\n\
	each object in list 2 with separation < d.\n\
	Used for frame registration.\n\
	Issues a warning if > 1 neighbour is found.\n\
	Input files must have x,y in 1st two columns.\n\
	Outputs\n\
		line1 line2 r_12 \n\
	where r_12 = separation and line1 = line read from 1st file.\n\
	Options are:\n\
		-d d	# maximum separation (3)\n\n"


#define MAX_STRING_LENGTH 1024

typedef struct particle {
        double                  x, y;
        struct  particle        *next;
	char	line[MAX_STRING_LENGTH];
} particle;

#define TINY	1.e-10
void	getcoords(particle *theparticle, int *i, int *j);
double	distance(particle *p1, particle *p2);

double		xmin, xmax, ymin, ymax, d;

main(int argc, char *argv[])
{
	int		arg, Nx, Ny, ix, iy, ix2, iy2, len;
	double		x, y, dist, distmin;
	particle	***grid, *baseparticle = NULL, *theparticle, *newparticle, *nextparticle;
	particle	*closest, second;
	FILE		*ipf;
	char		line[MAX_STRING_LENGTH];

	/* defaults */
	d = 3.0;

	/* parse args */
	if (argc < 3) {
		fprintf(stderr, usage);
		exit(-1);
	}
	arg = 3;
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			case 'd':
				sscanf(argv[arg++], "%lf", &d);
				break;
			default:
				fprintf(stderr, usage);
				exit(-1);
		}	
	}


	/* first check we can open 2nd list, and then read objects from 1st */
	/* into list headed by baseparticle and also figure xmin, xmax */
	if (!(ipf = fopen(argv[2], "r"))) {
		fprintf(stderr, "mergelists: cannot open %s\n", argv[2]);
		exit(-1);
	}
	fclose(ipf);
	if (!(ipf = fopen(argv[1], "r"))) {
		fprintf(stderr, "mergelists: cannot open %s\n", argv[1]);
		exit(-1);
	}
	while (fgets(line, 1024, ipf)) {
		if (line[0] == '#')
			continue;
		if (2 != sscanf(line, "%lf %lf", &x, &y)) {
			fprintf(stderr, "mergelists: format error in file %s\n", argv[1]);
			exit(-1);
		}
		newparticle = (particle *) calloc(1, sizeof(particle));
		newparticle->x = x;
		newparticle->y = y;
		len = strlen(line);
		line[len - 1] = '\0';
		strcpy(newparticle->line, line);
		newparticle->next = baseparticle;
		baseparticle = newparticle;
		if (x > xmax) xmax = x;
		if (x < xmin) xmin = x;
		if (y > ymax) ymax = y;
		if (y < ymin) ymin = y;
	}
	fclose(ipf);

	Nx = ceil ((xmax - xmin) / d);
	Ny = ceil ((ymax - ymin) / d);

	/* create the grid of particle ptrs */
	grid = (particle ***) calloc(Nx, sizeof(particle **));
	for (ix = 0; ix < Nx; ix++) {
		grid[ix] = (particle **) calloc(Ny, sizeof(particle *));
	}

	/* install the particles in the grid */
	theparticle = baseparticle;
	while (theparticle) {
		getcoords(theparticle, &ix, &iy);
		if (ix < 0 || ix >= Nx || iy < 0 || iy >= Ny) {
			theparticle = theparticle->next;
			continue;
		}
		nextparticle = theparticle->next;
		theparticle->next = grid[ix][iy];
		grid[ix][iy] = theparticle;
		theparticle = nextparticle;
	}

	/* now we process the 2nd list of particles */
	ipf = fopen(argv[2], "r");
	while (fgets(line, 1024, ipf)) {
		if (line[0] == '#')
			continue;
		if (2 != sscanf(line, "%lf %lf", &(second.x),  &(second.y))) {
			fprintf(stderr, "mergelists: format error in file %s\n", argv[2]);
			exit(-1);
		}
		len = strlen(line);
		line[len - 1] = '\0';
		strcpy(second.line, line);
		getcoords(&second, &ix2, &iy2);
		for (ix = ix2 - 1; ix <= ix2 + 1; ix++) {
			if (ix < 0 || ix >= Nx)
				continue;
			for (iy = iy2 - 1; iy <= iy2 + 1; iy++) {
				if (iy < 0 || iy >= Ny)
					continue;
				theparticle = grid[ix][iy];
				distmin = 10 * d;
				while (theparticle) {
					dist = distance(theparticle, &second);
					if (dist < d && distmin < d)
						fprintf(stderr, "mergelists: warning multiple neighbours\n");
					if (dist < d && dist < distmin) {
						distmin = dist;
						closest = theparticle;
					}
					theparticle = theparticle->next;
				}
				if (distmin < 10 * d) {
					fprintf(stdout, "%s %s %10.5lf\n", 
						closest->line, second.line, distmin);
				}
			}
		}
	}
}

void	getcoords(particle *theparticle, int *ix, int *iy)
{
	*ix = floor((theparticle->x - xmin) / d);
	*iy = floor((theparticle->y - ymin) / d);	
}



double	distance(particle *p1, particle *p2)
{
	return (sqrt((p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y)));
}

