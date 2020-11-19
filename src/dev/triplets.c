/*
 * triplets.c
 *
 */


#include <stdio.h>
#include <math.h>

#define usage "usage: triplets [dmax]\n\
	reads a list of coords and finds triangles with max side dmax\n\
	orders vertices so long, medium and short sides are\n\
	L = AB; M = BC; S = CA and prints out\n\
	L, M, S, phi, xA, yA where\n\
	phi = atan2(yB - yA, xB - xA)\n\n"

#define		WHOPPER		1.e10
#define 	MAX_OBJECTS	1000

#define	dist(i,j)	sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]))

double	x[MAX_OBJECTS], y[MAX_OBJECTS];
void 	print(int ia, int ib, int ic);

main(int argc, char *argv[])
{
	double	dmax = WHOPPER;
	char	line[1024];
	int	n, i0, i1, i2, ntriplets;
	double	d01, d12, d20;

	if (argc == 2 && argv[1][0] == '-') {
		fprintf(stderr, "%s", usage);
		exit(-1);
	}
	if (argc == 2)
		sscanf(argv[1], "%lf", &dmax);

	/* read the objects */
	n = 0;
	while (fgets(line, 1024, stdin)) {
		if (line[0] == '#')
			continue;
		sscanf(line, "%lf %lf ", x + n, y + n);
		n++;
	}

	/* find the triplets */
	fprintf(stdout, "#        L          M          S        phi         xA         yA\n");
	ntriplets = 0;
	for (i0 = 0; i0 < n; i0++) {
		for (i1 = i0 + 1; i1 < n; i1++) {
			d01 = dist(i0,i1);
			if (d01 > dmax)
				continue;
			for (i2 = i1 + 1; i2 < n; i2++) {
				d12= dist(i1,i2);
				if (d12 > dmax)
					continue;
				d20= dist(i2,i0);
				if (d20 > dmax)
					continue;
				if (d01 == d12 || d12 == d20 || d20 == d01)
					continue;
				if (inorder(d01,d12,d20)) {
					print(i0,i1,i2);
					continue;
				}
				if (inorder(d01,d20,d12)) {
					print(i1,i0,i2);
					continue;
				}
				if (inorder(d12,d01,d20)) {
					print(i2,i1,i0);
					continue;
				}
				if (inorder(d12,d20,d01)) {
					print(i1,i2,i0);
					continue;
				}
				if (inorder(d20,d01,d12)) {
					print(i2,i0,i1);
					continue;
				}
				if (inorder(d20,d12,d01)) {
					print(i0,i2,i1);
					continue;
				}
			}
		}
	}
}


void print(int ia, int ib, int ic)
{
	double	phi;
	
	phi = atan2(y[ib]-y[ia], x[ib]-x[ia]);
	fprintf(stdout, "%10.3lf %10.3lf %10.3lf %10.5lf %10.3lf %10.3lf\n",
		dist(ia,ib), dist(ib,ic), dist(ic,ia), phi, x[ia], y[ia]);
}


int inorder(double d1, double d2, double d3)
{
	return (d1 > d2 && d2 > d3);
}
