#define usage "\n\n\n\
NAME\n\
	unpackextensions --- extract extensions from multi-part FITS file\n\
\n\
SYNOPSIS\n\
	unpackextensions [-n nextensions] [-N NAXIS] extname fmt [base]\n\
\n\
DESCRIPTION\n\
	\"unpackextensions\" is used to separate a multi-component\n\
	FITS file where a number of separate images (of arbitrary\n\
	dimensionality) have been packed as extensions.\n\
\n\
	We first read the primary header, check that it has \n\
	EXTEND = T and get the number of extensions from the\n\
	NEXTEND header item value.\n\
	Then, for each extension, we read the header and image data and\n\
	write them to a file with name 'filename' given by\n\
		sprintf(filename, fmt, EXTNAME) where\n\
	'fmt' is a format string and EXTNAME is the integerized value\n\
	of the header item named 'extname'.\n\
\n\
	For example, to unpack a CFHT archive format CFH12K file\n\
	504338o.fits to 504338o00.fits 504338o01.fits ... 504338o11.fits\n\
	you can do\n\
\n\
		unpackextensions IMAGEID 504338o\\%.2d.fits < 504338o.fits\n\
\n\
	since IMAGEID (or CHIPID) contains the chip number as an integer.\n\
\n\
	If extname = index1D then the filename is generated as\n\
		sprintf(filename, fmt, i)\n\
	where i = 0, .... n-1 is an internally generated index.\n\
\n\
	If extname = index2D then the integer argument base must be supplied\n\
	and the filename is generated as\n\
		printf(filename, fmt, j, k)\n\
	where k = i % base and j = (i - k) / base.\n\
\n\
	if extname = index2d then the integer argument base must be supplied\n\
	and the filename is generated as\n\
		printf(filename, fmt, j, k)\n\
	where k = base - 1 - i % base and j = (i - k) / base.\n\
	This may be useful for unpacking GPC data.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@hawaii.edu\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/args.h"

int		main(int argc, char *argv[])	
{
	fitsheader	*fits;
	fitscomment	*thecomment;
	char		*extname, *fmt, *flag, filename[128];
	int		i, id, nbytes, dim, rem, nim, base, nextlim = 0, naxis3 = 0;
	FILE		*opf;
	void		*data;

	/* parse args */
	argsinit(argc, argv, usage);
	while (nextargtype() == FLAG_ARG) {
		flag = getflag();
		switch (flag[0]) {
			case 'n':
				nextlim = getargi();
				break;
			case 'N':
				naxis3 = getargi();
				break;
			default:
				error_exit(usage);
				break;
		}
	}
	extname = getargs();
	fmt = getargs();
	if (!strcmp(extname, "index2D") || !strcmp(extname, "index2d")) {
		base = getargi();
	}
	fits = readfitsheader(stdin);
	if (!(fits->hasextensions)) {
		error_exit("unpackextensions: primary header contains no extensions\n");
	}
	nim = (nextlim ? nextlim : fits->nextensions);

	for (i = 0; i < nim; i++) {
		fits = readfitsheader(stdin);
		if (!naxis3 || naxis3 == fits->ndim) {
			if (strncmp(extname, "index", 5)) {
				thecomment = getcommentbyname(extname, fits);
				if (thecomment) {
					id = (int) getnumericvalue(thecomment);
				} else {
					error_exit("unpackextensions: couldn't find extname header item\n");
				}
				sprintf(filename, fmt, id);
			} else {
				if (!strcmp(extname, "index1D")) {
					sprintf(filename, fmt, i);
				} else {
					if (!strcmp(extname, "index2D")) {
						sprintf(filename, fmt, (i - i % base) / base, i % base);
					} else {
						if (!strcmp(extname, "index2d")) {
							sprintf(filename, fmt, (i - i % base) / base, base - 1 - i % base);
						} else {
							error_exit("unpackextensions: can't understand extname argument\n");
						}
					}
				}
			}
			/* fprintf(stderr, "# filename = %s\n", filename); */
			opf = fopen(filename, "w");
			if (!opf) {
				error_exit("unpackextensions: couldn't open output file\n");
			}
			fits->isextension = 0;
			fits->opstream = opf;
			fits->opbyteorder = BIG_ENDIAN_BYTE_ORDER;
			writefitsheader(fits);
		}
		nbytes = pixsize(fits->extpixtype);
		for (dim = 0; dim < fits->ndim; dim++) {
			nbytes *= fits->n[dim];
		}
		rem = nbytes % 2880;
		if (rem) {
			nbytes += (2880 - rem);
		}
		data = calloc(nbytes, sizeof(char));
		fread(data, sizeof(char), nbytes, stdin);
		if (!naxis3 || naxis3 == fits->ndim) {
			fwrite(data, sizeof(char), nbytes, opf);
		}
		free(data);
	}

	exit(0);
}

