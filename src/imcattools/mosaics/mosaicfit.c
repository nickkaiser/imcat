/*
 * mosaicfit.c
 *
 */

#define usage "\n\
NAME\n\
	mosaicfit --- fit for transformation coefficients for a set of exposures\n\
\n\
SYNOPSIS\n\
	mosaicfit [options...]\n\
\n\
DESCRIPTION\n\
	'mosaicfit' reads a catalogue containing the result of\n\
	merging all pairs of 'mosaiccat's (as created by 'mergemosaiccats1'\n\
	or 'mergemosaiccats2') and which must contain entries for the\n\
	spatial coords 'x', chip-number 'chip' and exp number 'exp'\n\
	and magnitude 'mag' from stdin, and fits a linear\n\
	model in which x,y are related to \"detector coords\" xd,yd by\n\
\n\
		xd = x + phi_c * y + dx_c\t\n\
		yd = y - phi_c * x + dy_c\t\n\
\n\
	so chips are rotated thru phi_c and displaced by dx_c, dy_c\n\
	relative to coordinate frame define by chip-0\n\
	and where sky coords X_e,Y_e are related to detector coords by\n\
\n\
		X_e = (1 + alpha * (xd * xd + yd * yd)) xd\t\n\
		Y_e = (1 + alpha * (xd * xd + yd * yd)) yd\t\n\
\n\
	and where sky coords in frame defined by exposure-0 are\n\
\n\
		X = X_e + Phi00_e * X_e + Phi01_e * Y_e + dX_e\t\n\
		Y = Y_e + Phi10_e * X_e + Phi11_e * Y_e + dY_e\t\n\
\n\
	which allows for pointing shifts and rotations as well as\n\
	any scale change or differential refraction.\n\
	The model is linearised - so only valid for small\n\
	alpha, phi_c, dx_c, dy_c, Phi_e (dX_e, dY_e can be large though)\n\
	Solves by minimising squared residuals.\n\
	We also read magnitudes, which we model as:\n\
\n\
		m_ce = m + m_c + M_e\n\
\n\
	where m is the true magnitude and m_c and M_e are magnitude\n\
	offsets for chip and exposure (relative to chip-0, exp-0).\n\
	Outputs coefficients in tabular form:\n\
		alpha\t\n\
		0	0	0	0	0	0	0\t\n\
			....\t\n\
		Phi00_m	Phi01_m	Phi10_m	Phi11_m	dX_m	dY_m	M_m\t\n\
			....\t\n\
		0	0	0	0\t\n\
			....\t\n\
		phi_n	dx_n	dy_n	m_n\t\n\
			....\t\n\
OPTIONS\n\
	Options are\n\
		-c Nc\t# number of chips (7)\n\
		-e Ne\t# number of exposures (11)\n\
		-n\t# don't compute magnitude shifts   \n\
\n\
</pre><p>See also <a href=\"./mosaicfitting.ps\"> mosaicfitting.ps </a><pre>\n\
\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <math.h>

#define SCALE 4096.0

main(int argc, char *argv[])
{
	int	arg, m, n, *indx, mcoefft, ncoefft, nobjects, c, cp, e, ep;
	double	**A, *B, d, *Cx, *Cy, **Am, *Bm, *Cm;
	double	x, y, xp, yp, rr, rrp, mag, magp;
	double	*Phi[2][2], *dX, *dY, *phi, *dx, *dy, alpha, *dM, *dm;
	char	line[1024];
	int	Nc, Ne, dXbase, dYbase, phibase, dxbase, dybase, dMbase, dmbase;
	int	Phi00base, Phi01base, Phi10base, Phi11base, domagshifts;
	FILE	*lcpipe;


	/* defaults */
	Nc = 7;
	Ne = 11;
	domagshifts = 1;
	
	/* parse args */
	arg = 1;
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			fprintf(stderr, usage);
			exit(-1);
		}
		switch (argv[arg++][1]) {
			case 'c':
				sscanf(argv[arg++], "%d", &Nc);
				break;
			case 'e':
				sscanf(argv[arg++], "%d", &Ne);
				break;
			case 'n':
				domagshifts = 0;
				break;
			default:
				fprintf(stderr, usage);
				exit(-1);
				break;
		}
	}

	ncoefft = 1 + 6 * (Ne - 1) + 3 * (Nc - 1);
        B = (double *) calloc(ncoefft, sizeof(double));
        Cx = (double *) calloc(ncoefft, sizeof(double));
        Cy = (double *) calloc(ncoefft, sizeof(double));
        A = (double **) calloc(ncoefft, sizeof(double *));
        for (m = 0; m < ncoefft; m++) {
                A[m] = (double *) calloc(ncoefft, sizeof(double));
        }
        indx = (int *) calloc(ncoefft, sizeof(int));


	if (domagshifts) {
		mcoefft = (Ne - 1) + (Nc - 1);
		Bm = (double *) calloc(mcoefft, sizeof(double));
		Cm = (double *) calloc(mcoefft, sizeof(double));
		Am = (double **) calloc(mcoefft, sizeof(double));
		for (m = 0; m < mcoefft; m++) {
			Am[m] =  (double *) calloc(mcoefft, sizeof(double));
		}
	}

	Phi[0][0] = (double *) calloc(Ne, sizeof(double));
	Phi[0][1] = (double *) calloc(Ne, sizeof(double));
	Phi[1][0] = (double *) calloc(Ne, sizeof(double));
	Phi[1][1] = (double *) calloc(Ne, sizeof(double));
	dX = (double *) calloc(Ne, sizeof(double));
	dY = (double *) calloc(Ne, sizeof(double));
	phi = (double *) calloc(Nc, sizeof(double));
	dx = (double *) calloc(Nc, sizeof(double));
	dy = (double *) calloc(Nc, sizeof(double));
	dm = (double *) calloc(Nc, sizeof(double));
	dM = (double *) calloc(Ne, sizeof(double));
	dMbase = -1;
	dmbase = Ne - 2;
	Phi00base =	0;
	Phi01base =	(Ne - 1);
	Phi10base =	2 * (Ne - 1);
	Phi11base =	3 * (Ne - 1);
	dXbase =	4 * (Ne - 1);
	dYbase =	5 * (Ne - 1);
	phibase =	6 * (Ne - 1);
	dxbase = 	6 * (Ne - 1) + (Nc - 1);
	dybase = 	6 * (Ne - 1) + 2 * (Nc - 1);

	if (!(lcpipe = popen("lc -o x chip exp mag", "r"))) {
		fprintf(stderr, "mosaicfit: unable to open lc-pipe for input\n");
		exit(-1);
	}
	nobjects = 0;
	while (fgets(line, 1024, lcpipe)) {
		if (line[0] == '#')
			continue;
		sscanf(line, "%lf %lf %lf %lf %d %d %d %d %lf %lf", &x, &y, &xp, &yp, &c, &cp, &e, &ep, &mag, &magp);
		x /= SCALE;
		y /= SCALE;
		xp /= SCALE;
		yp /= SCALE;
		if (c < 0 || c > (Nc - 1) || cp < 0 || cp > (Nc - 1)) {
			fprintf(stderr, "mosaicfit: chip number out of allowed range\n");
			exit(-1);
		}
		if (e < 0 || e > (Ne - 1) || ep < 0 || ep > (Ne - 1)) {
			fprintf(stderr, "mosaicfit: exposure number out of allowed range\n");
			exit(-1);
		}
		nobjects++;
		/* calculate C-vectors */
		/* first the alpha term */
		rr = x * x + y * y;
		rrp = xp * xp + yp * yp;
		Cx[0] = rr * x - rrp * xp;
		Cy[0] = rr * y - rrp * yp;
		/* now the exposure terms */
		for (m = 1; m < Ne; m++) {
			Cx[m + Phi00base] = Cx[m + Phi01base] = Cx[m + Phi10base] = Cx[m + Phi11base] = 0.0;
			Cx[m + dXbase] = Cx[m + dYbase] = 0.0;
			Cy[m + Phi00base] = Cy[m + Phi01base] = Cy[m + Phi10base] = Cy[m + Phi11base] = 0.0;
			Cy[m + dXbase] = Cy[m + dYbase] = 0.0;
			if (e == m) {
/*
				Cx[m + Phibase] += y;
				Cy[m + Phibase] -= x;
*/
				Cx[m + Phi00base] += x;
				Cx[m + Phi01base] += y;
				Cy[m + Phi10base] += x;
				Cy[m + Phi11base] += y;
				Cx[m + dXbase] += 1.0;
				Cy[m + dYbase] += 1.0;
			}
			if (ep == m) {
/*
				Cx[m + Phibase] -= yp;
				Cy[m + Phibase] += xp;
*/
				Cx[m + Phi00base] -= x;
				Cx[m + Phi01base] -= y;
				Cy[m + Phi10base] -= x;
				Cy[m + Phi11base] -= y;
				Cx[m + dXbase] -= 1.0;
				Cy[m + dYbase] -= 1.0;
			}
		}
		/* now the chip terms */
		for (n = 1; n < Nc; n++) {
			Cx[n + phibase] = Cx[n + dxbase] = Cx[n + dybase] = 0.0;
			Cy[n + phibase] = Cy[n + dxbase] = Cy[n + dybase] = 0.0;
			if (c == n) {
				Cx[n + phibase] += y;
				Cy[n + phibase] -= x;
				Cx[n + dxbase] += 1.0;
				Cy[n + dybase] += 1.0;
			}
			if (cp == n) {
				Cx[n + phibase] -= yp;
				Cy[n + phibase] += xp;
				Cx[n + dxbase] -= 1.0;
				Cy[n + dybase] -= 1.0;
			}
		}
		/* accumulate A matrix, B-vector */
		for (m = 0; m < ncoefft; m++) {
			B[m] -= (x - xp) * Cx[m] + (y - yp) * Cy[m];
			for (n = 0; n < ncoefft; n++) {
				A[m][n] += Cx[m] * Cx[n] + Cy[m] * Cy[n];
			}
		}
		if (domagshifts) {
			/* now we do the magnitude terms */
			/* first the exposure terms...*/
			for (m = 1; m < Ne; m++) {
				Cm[m + dMbase] = 0.0;
				if (e == m) {
					Cm[m + dMbase] += 1.0;
				}
				if (ep == m) {
					Cm[m + dMbase] -= 1.0;
				}
			}
			/* now the chip terms */
			for (n = 1; n < Nc; n++) {
				Cm[n + dmbase] = 0.0;
				if (c == n) {
					Cm[n + dmbase] += 1.0;
				}
				if (cp == n) {
					Cm[n + dmbase] -= 1.0;
				}
			}
			/* accumulate Am matrix, Bm-vector */
			for (m = 0; m < mcoefft; m++) {
				Bm[m] += (mag - magp) * Cm[m];
				for (n = 0; n < mcoefft; n++) {
					Am[m][n] += Cm[m] * Cm[n];
				}
			}
		}
	}

	if (nobjects <= ncoefft) {
		fprintf(stderr, "mosaicfit: too few objects\n");
		exit(-1);
	}

	myludcmp(A, ncoefft, indx, &d);
	mylubksb(A, ncoefft, indx, B);

	alpha = B[0] / (SCALE * SCALE);
	for (e = 1; e < Ne; e++) {
		Phi[0][0][e] = B[e + Phi00base];
		Phi[0][1][e] = B[e + Phi01base];
		Phi[1][0][e] = B[e + Phi10base];
		Phi[1][1][e] = B[e + Phi11base];
		dX[e] = SCALE * B[e + dXbase];
		dY[e] = SCALE * B[e + dYbase];
	}
	for (c = 1; c < Nc; c++) {
		phi[c] = B[c + phibase];
		dx[c] = SCALE * B[c + dxbase];
		dy[c] = SCALE * B[c + dybase];
	}

	if (domagshifts) {
		myludcmp(Am, mcoefft, indx, &d);
		mylubksb(Am, mcoefft, indx, Bm);

		for (e = 1; e < Ne; e++) {
			dM[e] = Bm[e + dMbase];
		}
		for (c = 1; c < Nc; c++) {
			dm[c] = Bm[c + dmbase];
		}
	}

	fprintf(stdout, "%13.8lg\n", alpha);
	for (e = 0; e < Ne; e++) {
		fprintf(stdout, "%13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg %13.8lg\n", 
			Phi[0][0][e], Phi[0][1][e], Phi[1][0][e], Phi[1][1][e], dX[e], dY[e], dM[e]);
	}
	for (c = 0; c < Nc; c++) {
		fprintf(stdout, "%13.8lg %13.8lg %13.8lg %13.8lg\n", phi[c], dx[c], dy[c], dm[c]);
	}
}
