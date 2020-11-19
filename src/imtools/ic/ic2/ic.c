#define usage "\n\n\
NAME\n\
		ic -- image calculator\n\
SYNOPSIS\n\
		ic [options....] rpnexpr fitsfile... \n\
		where options are:\n\
			-u 		# print usage\n\
			-c N1 N2	# create an image of this size\n\
			-f		# force float format output\n\
\n\
DESCRIPTION\n\
		\"ic\" does simple arithmetic on one or more images\n\
		according to the reverse-polish notation expression 'rpnexpr'.\n\
		Images are referred to in rpnexpr as '%%1', '%%2'.... and\n\
		must all have the same size.  If the fitsfilename is given\n\
		as '-' then that file will be read from stdin.\n\
		If any of the input images is MAGIC (as defined 'magic.h')\n\
		then the result will be MAGIC also.\n\
		The functions supported include all of the standard C math\n\
		library plus the operators '+', '-', '*', '/', and 'mult' \n\
		is provided as a synonym for '*' to avoid potential problems\n\
		if you invoke ic from a script. There are the logical\n\
		operations '>', '<', '>=', '<=', '!=' and the negation operator\n\
		'!'. There is a random number generator 'rand' and there are\n\
		two functions 'xp', 'yp' to get the horixontal and vertical pixel\n\
		positions respectively, and two further functions 'x', 'y' which\n\
		return the position in units of the lenght of the side of\n\
		the image.  There is a constant MAGIC (defined in magic.h - and\n\
		currently set to -32768) which is a flag for bad data.\n\
		There is a function 'if' (a.k.a. '?') which mimics the C\n\
		language '(c ? t : f)' which returns 't' or 'f' respectively\n\
		depending on the truth or falseness of the condition 'c'.  The\n\
		rpn syntax for this expression is 't f c ?' in which '?' pops the\n\
		condition 'c' followed by 'f' and then 't' and pushes 't' or 'f'\n\
		as appropriate. The condition 'c' will of course most likely be\n\
		a compound logical expression.  There is also a function 'enter'\n\
		which duplicates the top value of the stack.\n\
		Use -c option (with no input images) to generate an image from\n\
		scratch.  Use -f option to force output to be floating point.\n\
		Otherwise the output will have same pixtype as input, or, with\n\
		-c option, will have short (2 byte) pixels.\n\
		Output goes to stdout.\n\
\n\
EXAMPLES\n\
		Subtract b.fits from a.fits:\n\
\n\
			ic '%%1 %%2 -' a.fits b.fits\n\
\n\
		Take sqrt of image to be read from stdin; output as float:\n\
\n\
			ic -f '%%1 sqrt' -\n\
\n\
		Create a 512 x 512 image with a linear horizontal ramp and\n\
		multiply by 10:\n\
\n\
			ic -c 512 512 'x 10 *'\n\
\n\
		Replace all pixels in a.fits with value < 10 with MAGIC:\n\
\n\
			ic '%%1 MAGIC %%1 10 > ?' a.fits\n\
\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


#include "../../imlib/fits.h"
#include "error.h"
#include "../../imtools/magic.h"
#include  "getop.h"
#include  "operators.h"
#include  "ic.h"

#define	MAX_OPS	1024

int	nim = 0;
static	int	iy, N1, N2;
int	ix;
float	**f;

int		main(int argc, char *argv[])	
{
	int		arg = 1, M1, M2, im, iop, nop = 0;
	int		comc, comc2, pixtype, pixtype2;
	char		*comv[MAX_COMMENTS], *comv2[MAX_COMMENTS], *rpnexpr;
	float		*fout, op2, ff;
	FILE		**fitsf;
	int		usingstdin = 0, ok, createimage, floatoutput;
	op		*theop, *oplist[MAX_OPS];

	/* defaults */
	createimage = 0;
	floatoutput = 0;

	/* parse the args */
	if (argc < 2) {
		error_exit(usage);
	}
	while (arg < argc) {
                if (argv[arg][0] != '-') {	/* it must be the rpnexpr */
			break;
		} else {
                        switch (argv[arg++][1]) {
				case 'c':
					createimage = 1;
					sscanf(argv[arg++], "%d", &N1);
					sscanf(argv[arg++], "%d", &N2);
					comc = 0;
					break;
				case 'f':
					floatoutput = 1;
					break;
                                case 'u':
				default:
					error_exit(usage);
					break;
                        }
		}
	}
	if (argc == arg) {
		error_exit(usage);
	}

	rpnexpr = argv[arg++];
	
	/* now read the image headers and check they are same format */
	nim = argc - arg;
	if (nim && createimage) {
		error_exit(usage);
	}
	fitsf = (FILE **) calloc(nim, sizeof(FILE *));
	f = (float **) calloc(nim, sizeof(float *));
	for (im = 0; im < nim; im++) {
		if (strcmp(argv[arg], "-")) {
			if (!(fitsf[im] = fopen(argv[arg++], "r"))) {
				error_exit("ic: unable to open file for input\n");
			}
		} else {
			if (usingstdin) {
				error_exit("ic: Can't read two files from stdin\n");
			}
			fitsf[im] = stdin;
			usingstdin = 1;
			arg++;
		}
		set_fits_ipf(fitsf[im]);
		if (im == 0) {
			read_fits_head(&N1, &N2, &comc, comv);
			get_input_pixtype(&pixtype);
		} else {
			read_fits_head(&M1, &M2, &comc2, comv2);
			get_input_pixtype(&pixtype2);
			if ((M1 != N1) || (M2 != N2)|| (pixtype2 != pixtype)) {
				error_exit("ic: images have non-identical size or pixtype\n");
			}
		}
	}

	/* now parse the rpnstring */
	while (theop = getop(&rpnexpr)) {
		oplist[nop++]= theop;
		if (nop == MAX_OPS) {
			error_exit("ic: too many operators\n");
		}
	}

	if (floatoutput) {
		set_output_pixtype(FLOAT_PIXTYPE);
	} else {
		if (!createimage)
			set_output_pixtype(pixtype);
	}

	add_comment(argc, argv, &comc, comv);
	write_fits_head(N1, N2, comc, comv);

	/* allocate space for the image lines */
	fout = (float *) calloc(N1, sizeof(float));
	for (im = 0; im < nim; im++) {
		f[im] = (float *) calloc(N1, sizeof(float));
	}

	/* now do the biz... */
	for (iy = 0; iy < N2; iy++) {
		for (im = 0; im < nim; im++) {
			set_fits_ipf(fitsf[im]);
			fread_fits_line(f[im], N1);
		}
		for (ix = 0; ix < N1; ix++) {
			ok = 1;
			for (iop = 0; iop < nop; iop++) {
				theop = oplist[iop];
				funcarray[theop->type](theop);
			}
			if (ok) {
				fout[ix] = (float) pop();
			} else {
				fout[ix] = MAGIC;
/*				resetstack();*/
			}
		}
		fwrite_fits_line(fout, N1);		
	}
	exit(0);
}




double	xp(void)
{
	return ((double) ix);
}

double	yp(void)
{
	return ((double) iy);
}

double	x(void)
{
	return ((double) ix / (double) N1);
}

double	y(void)
{
	return ((double) iy / (double) N2);
}

