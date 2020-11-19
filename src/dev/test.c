#include <stdio.h>
#include <math.h>

double	drand48();

main()
{
	int	i, nstars, i1, j1, i2, j2;
	double	X, Y, x, y, r, R, alpha;

	nstars = 20;
	alpha = 2.e-10;
	for (i = 0; i < nstars; i++) {
		x = 2048 * drand48();
		y = 4096 * drand48();
		r = sqrt(x * x + y * y);
		R = r * (1 + alpha * r * r);
		X = x * R / r;
		Y = y * R / r;
		fprintf(stdout, "%10d %10d ", (int) floor(4096 - Y), (int) floor(2048 - X));
		x += 300;
		y += 1200;
		r = sqrt(x * x + y * y);
		R = r * (1 + alpha * r * r);
		X = x * R / r;
		Y = y * R / r;
		fprintf(stdout, "%10d %10d\n", (int) floor(4096 - Y), (int) floor(2048 - X));
	}	
}
