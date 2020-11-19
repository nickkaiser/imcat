/*
 * matchobjects.c
 *
 */


#include <stdio.h>
#include <math.h>

#define usage "usage: matchobjects.c [options...] a.tab b.tab\n\
	reads a pair of tables and tries to find matching lines.\n\
	Options are:\n\
		-a col prec	# absolute comparison on column col\n\
		-r col prec	# relative comparison on column col\n\
		-o col format	# output column col\n\
		-d col format	# output difference in column col\n\
	Default behaviour is equivalent to \'matchobjects -a 1 0 -o 2 \"%10.3\"\'\n\
	which outputs col2 from a.tab followed by col2 from b.tab provided\n\
	col1 values agree exactly\n\
	U
	Use -r to specify relative comparsion 
	and (dx,dy) is the translation in b-coords\n\n"

typedef struct triangle {
	double	L, M, S, phi, x, y, m, s;
	struct	triangle	*next;
} triangle;

#define	PI	M_PI

main(int argc, char *argv[])
{
	double		mmax, smax, L, M, S, phi, x, y, m, s, R, dphi, dx, dy, xx, yy;
	FILE		*af, *bf;
	triangle	*at, *atbase;
	char		line[1024];

	mmax = smax = 0.01;

	switch (argc) {
		case 5:
			sscanf(argv[3], "%lf", &mmax);
			sscanf(argv[4], "%lf", &smax);
			break;
		case 3:
			break;
		default:
			fprintf(stderr, "%s", usage);
			exit(-1);
	}

	/* open the triplets files */
	af = fopen(argv[1], "r");
	bf = fopen(argv[2], "r");
	if (!af || !bf) {
		fprintf(stderr, "matchpairs: failed to open input file\n");
	}
	
	/* read the first list into a null terminated linked list with head abaset */
	at=NULL;
	while (fgets(line, 1024, af)) {
		if (line[0] == '#')
			continue;
		atbase = (triangle *) calloc(1, sizeof(triangle));
		atbase->next = at;
		at = atbase;
		sscanf(line, "%lf %lf %lf %lf %lf %lf", 
			&(at->L), &(at->M), &(at->S),
			&(at->phi), &(at->x), &(at->y));
		at->m = at->M / at->L;
		at->s = at->S / at->L;
	}
	fclose(af);
	
	/* now process the second list */
	fprintf(stdout, "#        R       dphi         dx         dy\n");
	while (fgets(line, 1024, bf)) {
		if (line[0] == '#')
			continue;
		sscanf(line, "%lf %lf %lf %lf %lf %lf", 
			&L, &M, &S, &phi, &x, &y);
		m = M / L;
		s = S / L;
		at = atbase;
		while (at) {
			if (fabs(at->m - m) < mmax && fabs(at->s - s) < smax) {
				R = L / at->L;
				dphi = phi - at->phi;
				while (dphi > PI)
					dphi -= 2 * PI;
				while (dphi < -PI)
					dphi += 2 * PI;
				/* rotate x,y */
				xx = at->x * cos(dphi) - at->y * sin(dphi);
				yy = at->y * cos(dphi) + at->x * sin(dphi);
				dx = x - R * xx;
				dy = y - R * yy;
				fprintf(stdout, "%10.5lf %10.5lf %10.3lf %10.3lf\n",
					R, dphi, dx, dy);
			}
			at = at->next;
		}
	}	
}
