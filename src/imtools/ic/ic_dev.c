#define usage "\n\n\
NAME\n\
		ic -- image calculator\n\
SYNOPSIS\n\
		ic [options....] rpnexpr fitsfile... \n\
		where options are:\n\
			-u 		# print usage\n\
			-c N1 N2	# create an image of this size\n\
			-p pixtype	# specify output pixtype\n\
			-s seed		# seed rand num generator\n\
			-m magicsub	# substitute magic value\n\
\n\
DESCRIPTION\n\
		\"ic\" does simple arithmetic on one or more images\n\
		according to the reverse-polish notation expression 'rpnexpr'.\n\
		Images are referred to in rpnexpr as '%%1', '%%2'.... and\n\
		must all have the same size.  If the fitsfilename is given\n\
		as '-' then that file will be read from stdin.\n\
		If any of the input images is MAGIC (as defined 'magic.h')\n\
		then the result will be MAGIC also (though with '-m' option\n\
		'ic' will output 'magicsub' in place of the usual SHRT_MIN)\n\
		The functions supported include all of the standard C math\n\
		library plus the operators '+', '-', '*', '/', and 'mult' \n\
		is provided as a synonym for '*' to avoid potential problems\n\
		if you invoke ic from a script. There are the logical\n\
		operations '>', '<', '>=', '<=', '!=', '==' and the negation operator\n\
		'!'. There is a random number generator 'rand' which generates\n\
		a uniform random number on the range 0.0-1.0 and 'grand' which\n\
		generates a zero-mean, unit variance normal variate. There are\n\
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
		scratch.\n\
		Use -p option to specify output pixtpye which can be one of\n\
			 16	# signed 2-byte integer\n\
			-32	# 4-byte floating point\n\
			 32	# 4-byte int\n\
		Otherwise the output will have same pixtype as input, or, with\n\
		-c option, will have pixtype -32.\n\
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
#include "../../utils/arrays.h"
#include  "getop.h"
#include  "ic.h"

#define	MAX_OPS	1024
#define STACK_DEPTH 100
double  st[STACK_DEPTH];
int     pos = -1;

int	nim = 0;
static	int	ix, iy, N1, N2;
float	**f;

double	drand48();

int		main(int argc, char *argv[])	
{
	int		arg, M1, M2, im, iop, nop = 0;
	int		comc, comc2, pixtype, pixtype2, outputpixtype, badpix;
	char		*comv[MAX_COMMENTS], *comv2[MAX_COMMENTS], *rpnexpr, rpnstring[128], *argvcopy[100];
	float		*fout, op2, ff, magicsub, ***fin, **ffout;
	FILE		**fitsf;
	int		usingstdin = 0, ok, createimage, floatoutput;
	op		*theop, *oplist[MAX_OPS];
	long		seed;

	/* defaults */
	createimage = 0;
	outputpixtype = 0;
	magicsub = MAGIC;

	/* copy the args */
	for (arg = 0; arg < argc; arg++)
		argvcopy[arg] = argv[arg];
	arg = 1;

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
				case 'p':
					sscanf(argv[arg++], "%d", &outputpixtype);
					if (outputpixtype != SHORT_PIXTYPE &&
						outputpixtype != FLOAT_PIXTYPE &&
						outputpixtype != INT_PIXTYPE) {
							error_exit("ic: illegal pixtype\n");
					}
					set_output_pixtype(outputpixtype);
					break;
				case 's':
					sscanf(argv[arg++], "%ld", &seed);
					srand48(seed);
					break;
				case 'm':
					sscanf(argv[arg++], "%f", &magicsub);
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

	if (!outputpixtype && createimage)
		set_output_pixtype(outputpixtype = FLOAT_PIXTYPE);

	rpnexpr = argv[arg];
	sprintf(rpnstring, "'%s'", rpnexpr);
	argvcopy[arg++] = rpnstring;
	
	/* now read the image headers and check they are same format */
	nim = argc - arg;
	if (nim && createimage) {
		error_exit(usage);
	}
	fitsf = (FILE **) calloc(nim, sizeof(FILE *));
	f = (float **) calloc(nim, sizeof(float *));
	fin = (float ***) calloc(nim, sizeof(float **));
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
			fread_fits(fin + im, &N1, &N2, &comc, comv);
			get_input_pixtype(&pixtype);
			if (!outputpixtype)
				set_output_pixtype(outputpixtype = pixtype);
		} else {
			fread_fits(fin + im, &M1, &M2, &comc2, comv2);
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


	add_comment(argc, argvcopy, &comc, comv);
/*	write_fits_head(N1, N2, comc, comv); */

	/* allocate space for the output image */
	allocFloatArray(&ffout, N1, N2);
 
	/* allocate space for the image lines */
/*	fout = (float *) calloc(N1, sizeof(float));
	for (im = 0; im < nim; im++) {
		f[im] = (float *) calloc(N1, sizeof(float));
	}*/

	/* now do the biz... */
	for (iy = 0; iy < N2; iy++) {
		for (im = 0; im < nim; im++) {
			/*set_fits_ipf(fitsf[im]);
			fread_fits_line(f[im], N1);*/
			f[im] = fin[im][iy];
		}
		fout = ffout[iy];
		for (ix = 0; ix < N1; ix++) {
			badpix = 0;
			for (iop = 0; iop < nop; iop++) {
				theop = oplist[iop];
				switch (theop->type) {
					case CONSTANT_TYPE:
						pos++;
						if (pos >= STACK_DEPTH) {
							fprintf(stderr, "ic: st full\n");
							exit(-1);
						}
						st[pos] = theop->data;
						/* if (st[pos] == MAGIC)
							badpix = 1;*/
						break;
					case IM_VALUE_TYPE:
						pos++;
						if (pos >= STACK_DEPTH) {
							fprintf(stderr, "ic: st full\n");
							exit(-1);
						}
						st[pos] = f[theop->opno][ix];
						if (st[pos] == MAGIC)
							badpix = 1;
						break;
					case NUM0_FUNC_TYPE:
						pos++;
						if (pos >= STACK_DEPTH) {
							fprintf(stderr, "ic: st full\n");
							exit(-1);
						}
						switch (theop->opno) {
							case RAND:
								st[pos] = drand48();
								break;
							case X:
								st[pos] = x();
								break;
							case Y:
								st[pos] = y();
								break;
							case XP:
								st[pos] = xp();
								break;
							case YP:
								st[pos] = yp();
								break;
							case IF:
							case IFQ:
								pos -= 3;
								if (pos < 0) {
									fprintf(stderr, "ic: stack empty\n");
									exit(-1);
								}
								st[pos] = (double) (st[pos+2] ? st[pos] : st[pos+1]);
								break;
							case ENTER:
								st[pos] = st[pos-1];
								break;
							case GRAND:
								st[pos] = gasdev();
								break;
							default:
								fprintf(stderr, "ic: bad opno\n");
								exit(-1);
								break;
						}
						break;
					case NUM1_FUNC_TYPE:
						switch (theop->opno) {
							case EXP:
								st[pos] = exp(st[pos]);
								break;
							case LOG:
								st[pos] = log(st[pos]);
								break;
							case LOG10:
								st[pos] = log10(st[pos]);
								break;
							case FABS:
								st[pos] = fabs(st[pos]);
								break;
							case NOT:
								st[pos] = (double)(!((int) st[pos]));
								break;
							case SQRT:
								st[pos] = sqrt(st[pos]);
								break;
							case FLOOR:
								st[pos] = floor(st[pos]);
								break;
							case CEIL:
								st[pos] = ceil(st[pos]);
								break;
							case SIN:
								st[pos] = sin(st[pos]);
								break;
							case ASIN:
								st[pos] = asin(st[pos]);
								break;
							case SINH:
								st[pos] = sinh(st[pos]);
								break;
							case COS:
								st[pos] = cos(st[pos]);
								break;
							case ACOS:
								st[pos] = acos(st[pos]);
								break;
							case COSH:
								st[pos] = cosh(st[pos]);
								break;
							case TAN:
								st[pos] = tan(st[pos]);
								break;
							case ATAN:
								st[pos] = atan(st[pos]);
								break;
							case TANH:
								st[pos] = tanh(st[pos]);
								break;
							default:
								fprintf(stderr, "ic: bad opno\n");
								exit(-1);
								break;
						}
						break;
					case NUM2_FUNC_TYPE:
						pos--;
						if (pos < 0) {
							fprintf(stderr, "ic: st empty\n");
							exit(-1);
						}
						switch (theop->opno) {
							case TIMES0:
							case TIMES1:
								st[pos] = st[pos] * st[pos+1];
								break;
							case PLUS:
								st[pos] = st[pos] + st[pos+1];
								break;
							case MINUS:
								st[pos] = st[pos] - st[pos+1];
								break;
							case DIVIDE:
								st[pos] = st[pos] / st[pos+1];
								break;
							case POW:
								st[pos] = pow(st[pos], st[pos+1]);
								break;
							case GT:
								st[pos] = (double) (st[pos] > st[pos+1]);
								break;
							case GE:
								st[pos] = (double) (st[pos] >= st[pos+1]);
								break;
							case LT:
								st[pos] = (double) (st[pos] < st[pos+1]);
								break;
							case LE:
								st[pos] = (double) (st[pos] <= st[pos+1]);
								break;
							case EQ:
								st[pos] = (double) (st[pos] == st[pos+1]);
								break;
							case NE:
								st[pos] = (double) (st[pos] != st[pos+1]);
								break;
							case FMOD:
								st[pos] = fmod(st[pos], st[pos+1]);
								break;
							case ATAN2:
								st[pos] = atan2(st[pos], st[pos+1]);
								break;
							default:
								fprintf(stderr, "ic: bad opno\n");
								exit(-1);
								break;
						}
						break;
					default:
						fprintf(stderr, "ic: bad op type\n");
						exit(-1);
						break;
				}
			}
			fout[ix] = (badpix ? magicsub : st[pos]);
			pos--;
		}
/*		fwrite_fits_line(fout, N1); */		
	}
	fwrite_fits(ffout, N1, N2, comc, comv);
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

double gasdev(void)
{
        static int iset=0;
        static double gset;
        double fac,r,v1,v2;
        double ran0();
 
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
