/*
 * lintransformfit.c
 *
 */

#define usage "usage: lintransformfit [options...]\n\
\t\t-r\t\t# output residuals\n\
\treads file containing x, y, x', y' from stdin and fits a linear\n\
\ttransformation model:\n\
\t\tx' = x0 + Pxx x + Pxy y\n\
\t\ty' = y0 + Pyx x + Pyy y\n\
\tby least squares.\n\
\tOutputs dx, dy, Pxx, Pxy, Pyx, Pyy to stdout by default.\n\
\tUse -r option to output x, y, x', y', dx, dy where dx = x' - model.\n"

#include <stdio.h>
#include <math.h>

main(int argc, char *argv[])
{
	int	doresiduals, arg, m, n, *indx, ncoefft, nobjects;
	double	**A1, **A2, *B;
	double	dx, dy, x, y, xp, yp, x0, y0, Pxx, Pxy, Pyx, Pyy, d;
	double	sum1, sumx, sumy, sumxx, sumxy, sumyy;
	double	sumxp, sumyp, sumxpx, sumxpy, sumypx, sumypy;
	char	line[1024];
	FILE	*tempf;

	/* parse args */
	doresiduals = 0;
	arg = 1;
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			case 'r':
				doresiduals = 1;
				break;
			default:
				fprintf(stderr, usage);
				exit(-1);
				break;
		}
	}

	ncoefft = 3;
        B = (double *) calloc(ncoefft, sizeof(double));
        A1 = (double **) calloc(ncoefft, sizeof(double *));
        A2 = (double **) calloc(ncoefft, sizeof(double *));
        for (m = 0; m < ncoefft; m++) {
                A1[m] = (double *) calloc(ncoefft, sizeof(double));
                A2[m] = (double *) calloc(ncoefft, sizeof(double));
        }
        indx = (int *) calloc(ncoefft, sizeof(int));

	sum1 = sumx = sumy = sumxx = sumxy = sumyy = 0.0;
	sumxp = sumyp = sumxpx = sumxpy = sumypx = sumypy = 0.0;
	nobjects = 0;
	if (doresiduals)
		tempf = fopen("lintransformfit.tmp", "w");
	while (fgets(line, 1024, stdin)) {
		if (line[0] == '#')
			continue;
		sscanf(line, "%lf %lf %lf %lf", &x, &y, &xp, &yp);
		if (doresiduals)
			fprintf(tempf, "%lf %lf %lf %lf\n", x, y, xp, yp);
		nobjects++;
		sum1 += 1.0;
		sumx += x;
		sumy += y;
		sumxx += x * x;
		sumxy += x * y;
		sumyy += y * y;
		sumxp += xp;
		sumyp += yp;
		sumxpx += xp * x;
		sumxpy += xp * y;
		sumypy += yp * y;
		sumypx += yp * x;
	}
	fclose(tempf);
	if (nobjects < 3)
		exit(-1);

	A1[0][0] = A2[0][0] = sum1;
	A1[1][1] = A2[1][1] = sumxx;
	A1[2][2] = A2[2][2] = sumyy;
	A1[0][1] = A1[1][0] = A2[0][1] = A2[1][0] = sumx;
	A1[0][2] = A1[2][0] = A2[0][2] = A2[2][0] = sumy;
	A1[1][2] = A1[2][1] = A2[1][2] = A2[2][1] = sumxy;

	B[0] = sumxp;
	B[1] = sumxpx;
	B[2] = sumxpy;

	myludcmp(A1, ncoefft, indx, &d);
	mylubksb(A1, ncoefft, indx, B);

	x0 = B[0];
	Pxx = B[1];
	Pxy = B[2];

	B[0] = sumyp;
	B[1] = sumypx;
	B[2] = sumypy;

	myludcmp(A2, ncoefft, indx, &d);
	mylubksb(A2, ncoefft, indx, B);

	y0 = B[0];
	Pyx = B[1];
	Pyy = B[2];

	if (doresiduals) {
		tempf = fopen("lintransformfit.tmp", "r");
		for (n = 0; n < nobjects; n++) {
			fscanf(tempf, "%lf %lf %lf %lf", &x, &y, &xp, &yp);
			dx = xp - x0 - Pxx * x - Pxy * y;
			dy = yp - y0 - Pyx * x - Pyy * y;
			fprintf(stdout, "%10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf\n",
				x, y, xp, yp, dx, dy);
				
		}
		fclose(tempf);
		system("rm lintransformfit.tmp");
	} else {
		fprintf(stdout, "%10.3lf %10.3lf %10.5lf %10.5lf %10.5lf %10.5lf\n",
			x0, y0, Pxx, Pxy, Pyx, Pyy);
	}
}
