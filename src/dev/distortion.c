/*
 * distortion.c
 *
 * fit for quadratic field distortion in 8k data
 */

#define usage "usage: distortion i0 j0\n\twhere (i0, j0) is nominal origin for chip\n"

#include <stdio.h>
#include <math.h>

main(int argc, char *argv[])
{
	int	i0, j0, m, i1, i2, j1, j2, *indx, ncoefft;
	double	Zx, Zy, X1, X2, Y1, Y2, **A, *B, RR1, RR2;
	double	x, y, alpha, d;
	double	sum1, sumZx, sumZy, sumZZ, sumX, sumY, sumZR;
	char	line[1024];

	if (argc != 3) {
		fprintf(stderr, usage);
		exit(-1);
	}
	sscanf(argv[1], "%d", &i0);
	sscanf(argv[2], "%d", &j0);

	ncoefft = 3;
        B = (double *) calloc(ncoefft, sizeof(double));
        A = (double **) calloc(ncoefft, sizeof(double *));
        for (m = 0; m < ncoefft; m++) {
                A[m] = (double *) calloc(ncoefft, sizeof(double));
        }
        indx = (int *) calloc(ncoefft, sizeof(int));

	sum1 = sumZx = sumZy = sumZZ = sumX = sumY = sumZR = 0.0;
	while (fgets(line, 1024, stdin)) {
		if (line[0] == '#')
			continue;
		sscanf(line, "%d %d %d %d", &i1, &j1, &i2, &j2);
		X1 = j1 + j0;
		Y1 = i1 + i0;
		X2 = j2 + j0;
		Y2 = i2 + i0;
		RR1 = (X1 * X1 + Y1 * Y1);
		RR2 = (X2 * X2 + Y2 * Y2);
		Zx = RR2 * X2 - RR1 * X1;
		Zy = RR2 * Y2 - RR1 * Y1;
		sum1	+= 1.0;
		sumZx	+= Zx;
		sumZy 	+= Zy;
		sumZZ	+= Zx * Zx + Zy * Zy;
		sumX	+= X2 - X1;
		sumY	+= Y2 - Y1;
		sumZR	+= Zx * (X2 - X1) + Zy * (Y2 - Y1);
	}

	A[0][0] = A[1][1] = sum1;
	A[0][2] = A[2][0] = sumZx;
	A[2][1] = A[1][2] = sumZy;
	A[2][2] = sumZZ;
	B[0] = sumX;
	B[1] = sumY;
	B[2] = sumZR;

	myludcmp(A, ncoefft, indx, &d);
	mylubksb(A, ncoefft, indx, B);

	x = B[0];
	y = B[1];
	alpha = B[2];

	fprintf(stdout, "x = %10.3lf\ny = %10.3lf\nalpha = %10g\n", x, y, alpha);
}
