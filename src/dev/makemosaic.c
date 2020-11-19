#include <stdio.h>
#include <math.h>

double	drand48();

#define MAXOBJS 10000
#define	NX	4
#define NY	2
#define N1	2048
#define N2	4096
#define L	10240
#define	GAP	50

main()
{
	int	frame, i, ntot, ix, iy, X, Y;
	double	x[MAXOBJS], y[MAXOBJS], dx[NX][NY], dy[NX][NY], phi[NX][NY], phimax, dmax;
	double	c, s, xc, yc, xx, yy, D;
	char	filename[128];
	FILE	*opf;

	ntot = 500;
	D = 1024;

	/* generate positions and rot angles */
	phimax = 0.02;
	dmax = 20;
	opf = fopen("mosaic_geometry.dat", "w");
	fprintf(opf, "#  X    Y        phi         dx         dy\n");
	for (X = 0; X < NX; X++) {
		for (Y = 0; Y < NY; Y++) {
			phi[X][Y] = 2 * phimax * (drand48() - 0.5);
			dx[X][Y] = X * (N1 + GAP) + 2 * dmax * (drand48() - 0.5);
			dy[X][Y] = Y * (N2 + GAP) + 2 * dmax * (drand48() - 0.5);
			fprintf(opf, "%4d %4d %10.5lf %10.3lf %10.3lf\n", X, Y, phi[X][Y], dx[X][Y], dy[X][Y]);
		}
	}
	fclose(opf);

	/* generate the objects in detector coords */
	for (i = 0; i < ntot; i++) {
		x[i] = L * drand48();
		y[i] = L * drand48();
	}

	/* now wewrite out the NX * NY lists */
	for (frame = 0; frame < 2; frame++) {
		for (X = 0; X < NX; X++) {
			for (Y = 0; Y < NY; Y++) {
				sprintf(filename, "mosaic%d_%d%d.list", frame, X, Y);
				opf = fopen(filename, "w");
				fprintf(opf, "#        x          y\n");
				c = cos(phi[X][Y]);
				s = sin(phi[X][Y]);
				for (i = 0; i < ntot; i++) {
					xx = x[i] - dx[X][Y] - D * frame;
					yy = y[i] - dy[X][Y] - D * frame;
					xc = c * xx + s * yy;
					yc = -s * xx + c * yy;
					if (0 < xc && xc < N1 && 0 < yc && yc < N2)
						fprintf(opf, "%10.3lf %10.3lf\n", xc, yc);
				}
				fclose(opf);
			}
		}
	}		
}
