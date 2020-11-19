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
			-h name val	# add header comment 'name = val'\n\
			-b BSCALE BZERO	# add scaling header information\n\
\n\
DESCRIPTION\n\
		\"ic\" does arithmetic on one or more images\n\
		according to the reverse-polish notation expression 'rpnexpr'.\n\
		Images are referred to in rpnexpr as '%%1', '%%2'.... and\n\
		must all have the same size.  If the fitsfilename is given\n\
		as '-' then that file will be read from stdin.\n\
		If fitsfilename is given as 'command |' then we read\n\
		from a pipe executing that command.\n\
\n\
		'ic' operates on images a line at a time, and so can be used\n\
		on very large images.  A reference to an input image such as\n\
		'%%1' causes a line of that image to be pushed onto a stack.\n\
		Single argument math functions operate on the line at the top\n\
		of the stack and multi-argument functions pop lines as\n\
		necessary and then push the resultant line.\n\
\n\
		If any of the input images is flagged as MAGIC (as defined 'fits.h')\n\
		then the result will be MAGIC also (though with '-m' option\n\
		'ic' will output 'magicsub' in place of the usual SHRT_MIN)\n\
\n\
		Fits header comments are inherited from the first image - other\n\
		comments are discarded.  You may add extra comments with the -h\n\
		option; see below.\n\
\n\
		The functions supported include all of the standard C math\n\
		library (including bessel functions j0(x), j1(x),\n\
		jn(n,x), y0(x), y1(x), yn(n,x))\n\
		plus the operators '+', '-', '*', '/', and 'mult' \n\
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
		a compound logical expression.  There are functions 'max' and 'min'\n\
		which pop two values and pushes the maximum or minimum respectively.\n\
\n\
		There is also a function 'enter' which duplicates the top value of\n\
		the stack, and a function 'swap' which unsurprisingly\n\
		swaps the top two values on the stack.\n\
\n\
		Use -c option (with no input images) to generate an image from\n\
		scratch.\n\
\n\
		Use -p option to specify output pixtpye which can be one of\n\
			  8	# 1-byte unsigned char\n\
			 16	# 2-byte signed integer\n\
			 32	# 4-byte signed int\n\
			-32	# 4-byte floating point\n\
			-64	# 8-byte floating point   \n\
		Otherwise the output will have same pixtype as that of the first\n\
		input image, or, with -c option, will have pixtype -32.\n\
\n\
		Use the -b option to apply scaling of pixel values on\n\
		output and record the BSCALE, BZERO values in the header.\n\
		Otherwise the BSCALE, BZERO values (if any) are inherited\n\
		from the (first) input image.\n\
		The definition of BSCALE and BZERO is such that the internal\n\
		values are computed from disk values as:\n\
\n\
			f_internal = BZERO + BSCALE * f_disk\n\
\n\
		Output goes to stdout.\n\
\n\
		The '-h' option can be used to add new header values or\n\
		to modify existing ones.  If an existing name is provided\n\
		then the existing value will be overwritten with the new value.\n\
		Otherwise, or if the name is 'HISTORY' or 'COMMENT', the new header line\n\
		will be appended.  The -h option can be given repeatedly to\n\
		insert a number of comments.\n\
\n\
EXAMPLES\n\
		Subtract b.fits from a.fits:\n\
\n\
			ic '%%1 %%2 -' a.fits b.fits\n\
\n\
		Take sqrt of image to be read from stdin; output as float:\n\
\n\
			ic -p -32 '%%1 sqrt' -\n\
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
		Filter to clip image at fmin = 0.25 fmax = 0.75:\n\
\n\
			ic '%%1 0.25 max 0.75 min' -\n\
\n\
AUTHOR\n\
		Nick Kaiser - kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


#include "../../imlib/fits.h"
#include "error.h"
#include "utils/iostream.h"
#include  "getop.h"
#include  "ic.h"

#define LOOP for(ix=0;ix<N1;ix++)

#define	MAX_OPS	1024
#define STACK_DEPTH 100
double  *st[STACK_DEPTH];
int     pos;

int	nim = 0;
static	int	ix, iy, N1, N2;
float	**f;

#define MAX_NEW_COMMENTS 100

#define	MAGIC FLOAT_MAGIC

double	drand48();

int		main(int argc, char *argv[])	
{
	int		arg, argccopy, M1, M2, im, iop, nop = 0, idim;
	int		pixtype, pixtype2, outputpixtype, *badpix;
	char		*rpnexpr, rpnstring[128], *argvcopy[100];
	float		*fout;
	float		op2, ff, magicsub;
	iostream	**ipstream;
	int		usingstdin = 0, ok, createimage, floatoutput;
	op		*theop, *oplist[MAX_OPS];
	long		seed;
	double		*tempdblptr;
	char		*newcommentname[MAX_NEW_COMMENTS], *newcommentval[MAX_NEW_COMMENTS];
	int		ncomments, icom, bscaling;
	double          bscale,  bzero;
	fitsheader	**fits, *fitsout;
	fitscomment	*com;

	/* defaults */
	createimage = 0;
	outputpixtype = 0;
	magicsub = MAGIC;
	ncomments = 0;
	bscaling = 0;

	/* start the args copy */
	arg = argccopy = 0;
	argvcopy[argccopy++] = argv[arg++];

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
 					argvcopy[argccopy++] = argv[arg-1];
					argvcopy[argccopy++] = argv[arg];
					sscanf(argv[arg++], "%d", &N1);
					argvcopy[argccopy++] = argv[arg];
					sscanf(argv[arg++], "%d", &N2);
					break;
				case 'p':
 					argvcopy[argccopy++] = argv[arg-1];
					argvcopy[argccopy++] = argv[arg];
					sscanf(argv[arg++], "%d", &outputpixtype);
					if (outputpixtype != UCHAR_PIXTYPE &&
						outputpixtype != SHORT_PIXTYPE &&
						outputpixtype != FLOAT_PIXTYPE &&
						outputpixtype != DBL_PIXTYPE &&
						outputpixtype != INT_PIXTYPE) {
							error_exit("ic: illegal pixtype\n");
					}
					break;
				case 's':
 					argvcopy[argccopy++] = argv[arg-1];
					argvcopy[argccopy++] = argv[arg];
					sscanf(argv[arg++], "%ld", &seed);
					srand48(seed);
					break;
				case 'm':
 					argvcopy[argccopy++] = argv[arg-1];
					argvcopy[argccopy++] = argv[arg];
					sscanf(argv[arg++], "%f", &magicsub);
					break;
				case 'h':
					newcommentname[ncomments] = argv[arg++];
					newcommentval[ncomments]  = argv[arg++];
					ncomments++;
					break;
				case 'b':
					bscaling = 1;
 					argvcopy[argccopy++] = argv[arg-1];
					argvcopy[argccopy++] = argv[arg];
					sscanf(argv[arg++], "%lf", &bscale);
					argvcopy[argccopy++] = argv[arg];
					sscanf(argv[arg++], "%lf", &bzero);
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
		outputpixtype = FLOAT_PIXTYPE;

	rpnexpr = argv[arg];
	sprintf(rpnstring, "'%.68s'", rpnexpr);
	argvcopy[argccopy++] = rpnstring;
	arg++;

	/* now read the image headers and check they are same size */
	nim = argc - arg;
	if (nim && createimage) {
		error_exit(usage);
	}
	ipstream = (iostream **) calloc(nim, sizeof(iostream *));
	fits = (fitsheader **) calloc(nim, sizeof(fitsheader *));
	f = (float **) calloc(nim, sizeof(float *));
	for (im = 0; im < nim; im++) {
		ipstream[im] = openiostream(argv[arg++], "r");
		/*
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
		*/
		fits[im] = readfitsheader(ipstream[im]->f);
		if (im == 0) {
			N1 = fits[im]->n[0];
			N2 = 1;
			for (idim = 1; idim < fits[im]->ndim; idim++) {
				N2 *= fits[im]->n[idim];
			}
			if (!outputpixtype) {
				outputpixtype = fits[im]->extpixtype;
			}
		} else {
			M1 = fits[im]->n[0];
			M2 = 1;
			for (idim = 1; idim < fits[im]->ndim; idim++) {
				M2 *= fits[im]->n[idim];
			}
			if ((M1 != N1) || (M2 != N2)) {
				error_exit("ic: images have non-identical size or pixtype\n");
			}
		}
	}
	/* set internal value of image width */
	set_getop_N1(N1);
	/* allocate space for badpix */
	badpix = (int *) calloc(N1, sizeof(int));
	/* allocate the data space */
	for (pos = 0; pos < STACK_DEPTH; pos++) {
		st[pos] = (double *) calloc(N1, sizeof(double));
	}

	/* allocate output image */
	if (nim) {
		fitsout = copyfitsheader(fits[0]);
		fitsout->extpixtype = outputpixtype;
	} else {
		fitsout = new2Dfitsheader(N1, N2, outputpixtype);
	}
	if (bscaling) {
		fitsout->bscaling = 1;
		fitsout->bscale = bscale;
		fitsout->bzero = bzero;
	}
	fout = (float *) calloc(N1, sizeof(float));

	/* now parse the rpnstring */
	while (theop = getop(&rpnexpr)) {
		oplist[nop++]= theop;
		if (nop == MAX_OPS) {
			error_exit("ic: too many operators\n");
		}
	}

	/* and add any other comments */
	for (icom = 0; icom < ncomments; icom++) {
		if (strcmp(newcommentname[icom], "HISTORY") && strcmp(newcommentname[icom], "COMMENT")) {
			removenamedcomments(newcommentname[icom], fitsout);
		}
		com = newtextcomment(newcommentname[icom], newcommentval[icom], NULL);
		appendcomment(com, fitsout);
	} 

	/* add the history comment */
	add_comment(argccopy, argvcopy, fitsout);

	writefitsheader(fitsout);
	/* allocate space for the image lines */
	for (im = 0; im < nim; im++) {
		f[im] = (float *) calloc(N1, sizeof(float));
	}

	/* now do the biz... */
	pos = -1;
	for (iy = 0; iy < N2; iy++) {
		for (im = 0; im < nim; im++) {
			readfitsline(f[im], fits[im]);
		}
			LOOP { badpix[ix] = 0; };
			for (iop = 0; iop < nop; iop++) {
				theop = oplist[iop];
				switch (theop->type) {
					case CONSTANT_TYPE:
						pos++;
						if (pos >= STACK_DEPTH) {
							fprintf(stderr, "ic: st full\n");
							exit(-1);
						}
						LOOP { st[pos][ix] = theop->data[ix]; };
						/* if (st[pos] == MAGIC)
							badpix = 1;*/
						break;
					case IM_VALUE_TYPE:
						pos++;
						if (pos >= STACK_DEPTH) {
							fprintf(stderr, "ic: st full\n");
							exit(-1);
						}
						LOOP {
							/* st[pos][ix] = f[theop->opno][iy][ix]; */
							st[pos][ix] = f[theop->opno][ix];
							if (st[pos][ix] == MAGIC)
								badpix[ix] = 1;
						}
						break;
					case NUM0_FUNC_TYPE:
						pos++;
						if (pos >= STACK_DEPTH) {
							fprintf(stderr, "ic: st full\n");
							exit(-1);
						}
						switch (theop->opno) {
							case RAND:
								LOOP { st[pos][ix] = drand48(); }
								break;
							case X:
								LOOP { st[pos][ix] = x(); }
								break;
							case Y:
								LOOP { st[pos][ix] = y(); }
								break;
							case XP:
								LOOP { st[pos][ix] = xp(); }
								break;
							case YP:
								LOOP { st[pos][ix] = yp(); }
								break;
							case IF:
							case IFQ:
								pos -= 3;
								if (pos < 0) {
									fprintf(stderr, "ic: stack empty\n");
									exit(-1);
								}
								LOOP { st[pos][ix] = (double) (st[pos+2][ix] ? st[pos][ix] : st[pos+1][ix]); }
								break;
							case ENTER:
								LOOP { st[pos][ix] = st[pos-1][ix]; }
								break;
							case GRAND:
								LOOP { st[pos][ix] = gasdev(); }
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
								LOOP { st[pos][ix] = exp(st[pos][ix]); }
								break;
							case LOG:
								LOOP { st[pos][ix] = log(st[pos][ix]); }
								break;
							case LOG10:
								LOOP { st[pos][ix] = log10(st[pos][ix]); }
								break;
							case FABS:
								LOOP { st[pos][ix] = fabs(st[pos][ix]); }
								break;
							case NOT:
								LOOP { st[pos][ix] = (double)(!((int) st[pos][ix])); }
								break;
							case SQRT:
								LOOP { st[pos][ix] = sqrt(st[pos][ix]); }
								break;
							case FLOOR:
								LOOP { st[pos][ix] = floor(st[pos][ix]); }
								break;
							case CEIL:
								LOOP { st[pos][ix] = ceil(st[pos][ix]); }
								break;
							case SIN:
								LOOP { st[pos][ix] = sin(st[pos][ix]); }
								break;
							case ASIN:
								LOOP { st[pos][ix] = asin(st[pos][ix]); }
								break;
							case SINH:
								LOOP { st[pos][ix] = sinh(st[pos][ix]); }
								break;
							case COS:
								LOOP { st[pos][ix] = cos(st[pos][ix]); }
								break;
							case ACOS:
								LOOP { st[pos][ix] = acos(st[pos][ix]); }
								break;
							case COSH:
								LOOP { st[pos][ix] = cosh(st[pos][ix]); }
								break;
							case TAN:
								LOOP { st[pos][ix] = tan(st[pos][ix]); }
								break;
							case ATAN:
								LOOP { st[pos][ix] = atan(st[pos][ix]); }
								break;
							case TANH:
								LOOP { st[pos][ix] = tanh(st[pos][ix]); }
								break;
							case BESSJ0:
								LOOP { st[pos][ix] = j0(st[pos][ix]); }
								break;
							case BESSJ1:
								LOOP { st[pos][ix] = j1(st[pos][ix]); }
								break;
							case BESSY0:
								LOOP { st[pos][ix] = y0(st[pos][ix]); }
								break;
							case BESSY1:
								LOOP { st[pos][ix] = y1(st[pos][ix]); }
								break;
							case SWAP:
								if (pos < 1) {
									fprintf(stderr, "ic: stack empty\n");
									exit(-1);
								}
								tempdblptr = st[pos];
								st[pos] = st[pos-1];
								st[pos-1] = tempdblptr;
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
								LOOP { st[pos][ix] = st[pos][ix] * st[pos+1][ix]; }
								break;
							case PLUS:
								LOOP { st[pos][ix] = st[pos][ix] + st[pos+1][ix]; }
								break;
							case MINUS:
								LOOP { st[pos][ix] = st[pos][ix] - st[pos+1][ix]; }
								break;
							case DIVIDE:
								LOOP { st[pos][ix] = st[pos][ix] / st[pos+1][ix]; }
								break;
							case POW:
								LOOP { st[pos][ix] = pow(st[pos][ix], st[pos+1][ix]); }
								break;
							case GT:
								LOOP { st[pos][ix] = (double) (st[pos][ix] > st[pos+1][ix]); }
								break;
							case GE:
								LOOP { st[pos][ix] = (double) (st[pos][ix] >= st[pos+1][ix]); }
								break;
							case LT:
								LOOP { st[pos][ix] = (double) (st[pos][ix] < st[pos+1][ix]); }
								break;
							case LE:
								LOOP { st[pos][ix] = (double) (st[pos][ix] <= st[pos+1][ix]); }
								break;
							case EQ:
								LOOP { st[pos][ix] = (double) (st[pos][ix] == st[pos+1][ix]); }
								break;
							case NE:
								LOOP { st[pos][ix] = (double) (st[pos][ix] != st[pos+1][ix]); }
								break;
							case FMOD:
								LOOP { st[pos][ix] = fmod(st[pos][ix], st[pos+1][ix]); }
								break;
							case ATAN2:
								LOOP { st[pos][ix] = atan2(st[pos][ix], st[pos+1][ix]); }
								break;
							case BESSJN:
								LOOP { st[pos][ix] = jn((int) st[pos][ix], st[pos+1][ix]); }
								break;
							case BESSYN:
								LOOP { st[pos][ix] = yn((int) st[pos][ix], st[pos+1][ix]); }
								break;
							case MAX:
								LOOP { st[pos][ix] = (st[pos][ix] > st[pos+1][ix] ? st[pos][ix] : st[pos+1][ix]); }
								break;
							case MIN:
								LOOP { st[pos][ix] = (st[pos][ix] < st[pos+1][ix] ? st[pos][ix] : st[pos+1][ix]); }
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
			LOOP { fout[ix] = (badpix[ix] ? magicsub : st[pos][ix]); }
			writefitsline(fout, fitsout);
			pos--;
	}
	writefitstail(fitsout);
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
