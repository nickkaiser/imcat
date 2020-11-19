#define usage "\n\n\n\
NAME\n\
	deproject - deproject an image of an assumed spherically symmetric object\n\
\n\
SYNOPSIS\n\
	deproject rmin rmax nbins [options....] \n\
		-c xc yc	# spatial origin (N1/2 + 0.5, N2/2 + 0.5)\n\
		-s d		# spacing of points for volume integration (0.1)\n\
\n\
DESCRIPTION\n\
	\"deproject\" computes deprojection of an assumed spherically\n\
	symmetric structure from a FITS image in the style of\n\
	Fabian et al.\n\
\n\
	The model is\n\
		f_2D = int dz f_3D(sqrt(rp^2 + z^2)).\n\
\n\
	Input image is read from stdin through profile (which generates\n\
	an lc-format version of the asimuthal sum of the image F(r)\n\
	in nbins log-spaced shells ranging from rmin to rmax.\n\
\n\
	We then work inward from outermost shell, computing first the\n\
	volumes\n\
		V[ir] = 2 pi int int dz drp rp\n\
	where integrals are over volumes such that r = sqrt(rp^2 + z^2)\n\
	falls in some given bin r, and where we use spacing of d times r\n\
	to perform the integrals.\n\
\n\
	We then compute:\n\
		f_3D[i] = (F[i] - sum_j>i V[j] f[j]) / V[i]\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>

#include "../../utils/error.h"
#include "../../utils/ipbuff.h"



int		main(int argc, char *argv[])	
{
	int		arg = 4;
	double		rmin, rmax, d, r, dr, dlnr;
	FILE		*ipf, *opf;
	char		lcstring[128], astring[128];
	double		**F, FF, *V, z, rp, *f3D, *r1, *r2;
	int		nF, i, j, nbins;

	/* defaults */
	d = 0.01;

	if (argc < 4) {
		error_exit(usage);
	}
	if (1 != sscanf(argv[1], "%lf", &rmin))
		error_exit(usage);
	if (1 != sscanf(argv[2], "%lf", &rmax))
		error_exit(usage);
	if (1 != sscanf(argv[3], "%d", &nbins))
		error_exit(usage);

	dlnr = log(rmax / rmin) / nbins;
	
   	sprintf(lcstring, "profile -l %lf %lf %d ", rmin, rmax, nbins);

    	while (arg < argc) {
                if (*argv[arg] != '-')
                        error_exit(usage);
                switch (*(argv[arg++]+1)) {
                        case 'c':
				strcat(lcstring, "-c ");
				strcat(lcstring, argv[arg++]);
				strcat(lcstring, " ");
				strcat(lcstring, argv[arg++]);
 				strcat(lcstring, " ");
                               	break;
                        case 's':
				if (1 != sscanf(argv[arg++], "%lf", &d))
					error_exit(usage);
                                break;
                        default:
                                error_exit(usage);
                                break;
                }
        }
        
	/* read the summed profile */
	strcat(lcstring, "| lc -b -o fsum");

	ipf = popen(lcstring, "r");
	if (!ipf) {
		error_exit("deproject: failed to open input pipe\n");
	}
	F = readdoublebuff(1, ipf, &nF);
	if (nF != nbins) {
		error_exit("deproject: mismatch between nbins and number of records read\n");
	}

	V = (double *) calloc(nbins, sizeof(double));
	f3D = (double *) calloc(nbins, sizeof(double));
	r1 = (double *) calloc(nbins, sizeof(double));
	r2 = (double *) calloc(nbins, sizeof(double));
	for (i = 0; i < nbins; i++) {
		r1[i] = rmin * exp(i * dlnr);
		r2[i] = r1[i] * exp(dlnr);
	}

	opf = popen("lc -C -n r1 -n r2 -n f", "w");
	if (!opf) {
		error_exit("deproject: failed to open output pipe\n");
	}

	for (i = nbins - 1; i >= 0; i--) {
		/* compute volumes */
		for (j = 0; j < nbins; j++) {
			V[j] = 0;
		}
		dr = d * r1[i];
		for (rp = r1[i] + 0.5 * dr; rp < r2[i]; rp += dr) {
			for (z = 0.5 * dr; z < rmax; z += dr) {
				r = sqrt(rp * rp + z * z);
				j = (int) floor(log(r / rmin) / dlnr);
				if (j >= 0 && j < nbins) {
					V[j] += rp;
				}
			}
		}
		for (j = 0; j < nbins; j++) {
			V[j] *= 4 * M_PI * dr * dr;
		}
		FF = F[i][0];
		for (j = i + 1; j < nbins; j++) {
			FF -= V[j] * f3D[j];
		}
		f3D[i] = FF / V[i];
	}

	for (i = 0; i < nbins; i++) {
		fprintf(opf, "%13.8lg %13.8lg %13.8lg\n", r1[i], r2[i], f3D[i]);
	}
	pclose(opf);	
	exit(0);
}




