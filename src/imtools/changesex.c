#define usage "\n\n\n\
NAME\n\
	changesex - switch byte order of a FITS file\n\
\n\
SYNOPSIS\n\
	changesex [-u] [-n]\n\
\n\
DESCRIPTION\n\
	The filter \"changesex\" swaps the byte order in a fits image.\n\
\n\
	Imcat commands can work with either big- or little-endian images.\n\
	Since 2/99 imcat written images contain a header value\n\
	with keyword BYTEORDR whose value can be either BIG_ENDIAN or\n\
	LITTLE_ENDIAN, and the value of which is set by the presence\n\
	or absence of the environment variable IMCATSWAPFITSBYTES.\n\
	If IMCATSWAPFITSBYTES is set then the image is stored in\n\
	non-native byte order.\n\
\n\
	Non-imcat images are assumed to be BIG_ENDIAN, while pre 2/99\n\
	imcat images are assumed to be in the native byte order.\n\
	If you have set IMCATSWAPFITSBYTES then you will need to filter\n\
	old style images through 'changesex' to bring them into non-native order.\n\
\n\
	Use 'changesex -n' to force images with BYTEORDR=LITTLE_ENDIAN\n\
	to BIG_ENDIAN format for compatibility with non-imcat software.\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/args.h"
#include "changesex.h"

int		main(int argc, char *argv[])	
{
	int		N1, N2, nlines, dim, i;
	fitsheader	*fits;
	void		*f;
	char		*flag;
	int             forcebigendian = 0;

	/* parse args */
	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch(flag[0]) {
		case 'n':
		  forcebigendian = 1;
		  break;
		case 'u':
		default:
		  error_exit(usage); 
		}
	}

	fits = readfitsheader(stdin);
	N1 = fits->n[0];
	nlines = 1;
	for (dim = 1; dim < fits->ndim; dim++) {
		nlines *= fits->n[dim];
	}
	add_comment(argc, argv, fits);
	/* here's the problem.... we dont flip fits->opbyteorder, so
	   the writefitsheader routine assumes the default (what's there
	   already) even though it's actually doing the swapping.
	   just an oversight I guess */
	if(fits->ipbyteorder == BIG_ENDIAN_BYTE_ORDER && forcebigendian != 1){
	  fits->opbyteorder = LITTLE_ENDIAN_BYTE_ORDER;
	}
	else{
	  fits->opbyteorder = BIG_ENDIAN_BYTE_ORDER;
	}
	writefitsheader(fits);

	f = calloc(N1, pixsize(fits->extpixtype));
	for (i = 0; i < nlines; i++) {
	  fread(f, pixsize(fits->extpixtype), N1, stdin);
	  if(!(fits->ipbyteorder == BIG_ENDIAN_BYTE_ORDER && forcebigendian == 1)){ 
	    byteswapline(f, N1, pixsize(fits->extpixtype));
	  }

	    fwrite(f, pixsize(fits->extpixtype), N1, stdout);
	  
	}
	writefitstail(fits);
	exit(0);
}



