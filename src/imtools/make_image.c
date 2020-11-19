#define usage "\n\n\n\
NAME\n\
	make_image --- generates a mock CCD image\n\
\n\
SYNOPSIS\n\
	make_image ncolours N1 N2 [option...]\n\
		-o n r f0		# number, size and central SB of gaussian objects\n\
		-e n r f0		# number, size and central SB of exponential objects\n\
		-s sigma		# rms sky fluctuation\n\
		-b col			# a bad column\n\
		-r seed			# for ran num generator\n\
		-f			# generate float format image\n\
		-i			# generate 32 bit int format image\n\
\n\
DESCRIPTION\n\
	\"make_image\" generates a mock CCD image containing families\n\
	of gaussian profile objects + noise + bad column defects.\n\
	The first argument must be the number of colours.\n\
	The second argument must be the size of the image on a side.\n\
	Random number seed should precede -o or -s args\n\
	Cosmic rays are indicated by a negative gaussian scale length.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "make_image.h"
#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/arrays.h"

#define GAU	0
#define EXP 1

#define	MAGIC FLOAT_MAGIC

static	float 	**f[MAX_FITS_DIM];
static	int	ncolours;
static	int	seed;

double	drand48();

int		main(int argc, char *argv[])	
{
	int	Naxes, comc, N1, N2, n, arg = 1, bad_col = -1, colour;
	float	r, f0, sigma;
	fitsheader	*fits;
	fitscomment	*com;
	char	argstring[COM_LENGTH], *color[1];
	int	i, j;
	
	if (!check_args(argc, argv)) {
		error_exit(usage);
	}
	
	fits = new2Dfitsheader(0, 0, SHORT_PIXTYPE);

	if (1 != sscanf(argv[arg++], "%d", &ncolours)) {
		fprintf(stderr, "bad image size argument\n");
		error_exit(usage);
	}
	if (ncolours > 1) {
		fits->ndim = 3;
		fits->n[2] = ncolours;
	} else {
		fits->ndim = 2;
	}

	if (1 != sscanf(argv[arg++], "%d", &N1)) {
		fprintf(stderr, "bad image size argument\n");
		error_exit(usage);
	}
	fits->n[0] = N1;

	if (1 != sscanf(argv[arg++], "%d", &N2)) {
		fprintf(stderr, "bad image size argument\n");
		error_exit(usage);
	}
	fits->n[1] = N2;

	for (colour = 0; colour < ncolours; colour++) {
		allocFloatArray(&(f[colour]), N1, N2);
	}
		
	while (arg < argc) {
		if (*argv[arg] != '-') {
			fprintf(stderr, "expecting -opt\n");
			error_exit(usage);
		}
		switch (*(argv[arg++]+1)) {
			case 'o':
				if (1 != sscanf(argv[arg++], "%d", &n)) {
					fprintf(stderr, "bad object number argument\n");
					error_exit(usage);
				}
				if (1 != sscanf(argv[arg++], "%f", &r)) {
					fprintf(stderr, "bad object size argument\n");
					error_exit(usage);
				}
				if (1 != sscanf(argv[arg++], "%f", &f0)) {
					fprintf(stderr, "bad object brightness argument\n");
					error_exit(usage);
				}
				add_population(N1, N2, n, f0, r, GAU);
				break;
			case 'e':
				if (1 != sscanf(argv[arg++], "%d", &n)) {
					fprintf(stderr, "bad object number argument\n");
					error_exit(usage);
				}
				if (1 != sscanf(argv[arg++], "%f", &r)) {
					fprintf(stderr, "bad object size argument\n");
					error_exit(usage);
				}
				if (1 != sscanf(argv[arg++], "%f", &f0)) {
					fprintf(stderr, "bad object brightness argument\n");
					error_exit(usage);
				}
				add_population(N1, N2, n, f0, r, EXP);
				break;
			case 's':
				if (1 != sscanf(argv[arg++], "%f", &sigma)) {
					fprintf(stderr, "bad sky sigma argument\n");
					error_exit(usage);
				}
				add_noise(N1, N2, sigma);
				break;
			case 'b':
				if (1 != sscanf(argv[arg++], "%d", &bad_col)) {
					fprintf(stderr, "bad bad-column number argument\n");
					error_exit(usage);
				}
				add_bad_col(N1, N2, bad_col);
				break;
			case 'r':
				if (1 != sscanf(argv[arg++], "%d", &seed)) {
					fprintf(stderr, "bad seed argument\n");
					error_exit(usage);
				}
				srand48((long) seed);
				break;
			case 'f':
				fits->extpixtype = FLOAT_PIXTYPE;
				break;
			case 'i':
				fits->extpixtype = INT_PIXTYPE;
				break;
			default:
				return (0);
				break;
		}
	}
	
	add_comment(argc, argv, fits);

        com = newtextcomment("HISTORY", "", "");
	switch (ncolours) {
		case 1:
			sprintf(com->value, "color: gray");
			break;
		case 2:
			sprintf(com->value, "color: red green");
			break;
		case 3:
			sprintf(com->value, "color: red green blue");
			break;
		default:
			break;
	}
        appendcomment(com, fits);
        writefitsheader(fits);

	for (colour = 0; colour < ncolours; colour++) {
		for (i = 0; i < N2; i++) {
			writefitsline(f[colour][i], fits);
		}
	}
	writefitstail(fits);
	return (0);
}


int		check_args(int argc, char *argv[])
{
	int 	arg;
	
	if (argc < 4) {
		return (0);
	}
	arg = 3;
	if (*argv[arg++] == '-') {
		return (0);
	}
	while (arg < argc) {
		if (*argv[arg] != '-') {
			return (0);
		}
		switch (*(argv[arg]+1)) {
			case 'o':
				arg += 4;
				break;
			case 'e':
				arg += 4;
				break;
			case 's':
				arg += 2;
				break;
			case 'b':
				arg += 2;
				break;
			case 'r':
				arg += 2;
				break;
			case 'f':
				arg += 1;
				break;
			case 'i':
				arg += 1;
				break;
			default:
				return (0);
				break;
		}
		if (arg > argc) {
			return (0);
		}
	}
	return (1);
}





	
void	add_bad_col(int N1, int N2, int col)
{
	int j, colour;
	
	if (col < 0 || col >= N1)
		return;
	for (colour = 0; colour <ncolours; colour++)
		for (j = 0; j < N2; j++)
			f[colour][j][col] = MAGIC;
}





void	add_population(int N1, int N2, int n_objects, float csb, float r_obj, int type)
{
	float 	rr_obj, rr;
	int			i, j, ij_max, n, idum, ii, jj, i_obj, j_obj, colour;
	float	rgb[3];

	if (r_obj >= 0) {
		ij_max = 3 * r_obj;							/* only go out to 3-sigma */
		rr_obj = r_obj * r_obj;
		
		for (n = 0; n < n_objects; n++)
			{
				for (colour = 0; colour < ncolours; colour++)
					rgb[colour] = drand48();
				i_obj = floor(N2 * drand48());
				j_obj = floor(N1 * drand48());
				for (i = i_obj - ij_max; i <= i_obj + ij_max; i++)
					{
						ii = i;
						while (ii < 0)						/* for periodic bc's */
							ii += N2;
						while (ii >= N2)
							ii -= N2;
						for (j = j_obj - ij_max; j <= j_obj + ij_max; j++)
							{
								jj = j;
								while (jj < 0)
									jj += N1;
								while (jj >= N1)
									jj -= N1;
								rr = (i - i_obj) * (i - i_obj) + (j - j_obj) * (j - j_obj);
								for (colour = 0; colour < ncolours; colour++) {
									switch (type) {
										case GAU:
											f[colour][ii][jj] += rgb[colour] * csb * exp( - 0.5 * rr / rr_obj);
											break;
										case EXP:
											f[colour][ii][jj] += rgb[colour] * csb * exp( - sqrt(rr / rr_obj));
											break;
										default:
											error_exit("add_population: bad type");
											break;
									}
								}
							}
					}
			}
	}
	else
		for (n = 0; n < n_objects; n++){
			i = floor(N2 * drand48());
			j = floor(N1 * drand48());
			f[0][i][j] += csb;
		}
}





void		add_noise(int N1, int N2, float sigma)
{
	int 		i, j, idum, colour;

	for (colour = 0; colour < ncolours; colour++)
		for (i = 0; i < N2; i++)						
			for (j = 0; j < N1; j++)
				f[colour][i][j] += sigma * gasdev();
}





float gasdev(void)
{
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;
	float ran0();

	if  (iset == 0) {
		do {
			v1=2.0*drand48()-1.0;
			v2=2.0*drand48()-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}






