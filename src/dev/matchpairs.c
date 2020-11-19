/*
 * matchobjects.c
 *
 */


#include <stdio.h>
#include <math.h>

#define usage "usage: matchpairs a.prs b.prs [dprec [mprec]]\n\
	reads a pair of pair-files generated by \"pairs\" and finds\n\
	pairs of pairs whose separations (dr) match to precision dprec (1.0)\n\
	and whose mags and mag differences match to precision mprec (0.2).\n\
	We output dphi = phi_a - phi_b, dx, dy where latter are\n\
	translation in b-coord system.\n\n"

typedef struct pair {
	double 	x, y, m, dr, phi, dm;
	struct	pair	*next;
} pair;

#define	PI	M_PI

main(int argc, char *argv[])
{
	double		dprec, mprec;
	FILE		*af, *bf;
	pair		*pr, *prbase;
	char		line[1024];
	double 		x, y, m, dr, phi, dm, dphi, dx, dy, xx, yy;

	/* defaults */
	dprec = 1.0;
	mprec = 0.2;

	switch (argc) {
		case 5:
			sscanf(argv[4], "%lf", &mprec);
		case 4:
			sscanf(argv[3], "%lf", &dprec);
		case 3:
			break;
		default:
			fprintf(stderr, "%s", usage);
			exit(-1);
	}

	/* open the pairs files */
	af = fopen(argv[1], "r");
	bf = fopen(argv[2], "r");
	if (!af || !bf) {
		fprintf(stderr, "matchpairs: failed to open input file\n");
	}
	
	/* read the first list into a null terminated linked list with head prbase */
	pr=NULL;
	while (fgets(line, 1024, af)) {
		if (line[0] == '#')
			continue;
		prbase = (pair *) calloc(1, sizeof(pair));
		prbase->next = pr;
		pr = prbase;
		sscanf(line, "%lf %lf %lf %lf %lf %lf", 
			&(pr->x), &(pr->y), &(pr->m),
			&(pr->dr), &(pr->phi), &(pr->dm));
	}
	fclose(af);
	
	/* now process the second list */
	fprintf(stdout, "#     dphi         dy         dy\n");
	while (fgets(line, 1024, bf)) {
		if (line[0] == '#')
			continue;
		sscanf(line, "%lf %lf %lf %lf %lf %lf", 
			&x, &y, &m, &dr, &phi, &dm);
		pr = prbase;
		while (pr) {
			if (fabs(pr->dr - dr) < dprec && fabs(pr->m - m) < mprec && fabs(pr->dm - dm) < mprec) {
				dphi = phi - pr->phi;
				while (dphi > PI)
					dphi -= 2 * PI;
				while (dphi < -PI)
					dphi += 2 * PI;
				/* rotate x,y */
				xx = pr->x * cos(dphi) - pr->y * sin(dphi);
				yy = pr->y * cos(dphi) + pr->x * sin(dphi);
				dx = x - xx;
				dy = y - yy;
				fprintf(stdout, "%10.5lf %10.3lf %10.3lf\n",
					dphi, dx, dy);
			}
			pr = pr->next;
		}
	}	
}