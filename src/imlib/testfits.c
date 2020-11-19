#include <stdio.h>
#include <limits.h>
#include <math.h>
#include "fits.h"

main (int argc, char *argv[]) 
{
	fitsheader	*fits, *fitso;
	fitscomment	*com, *lastcom;
	int		dim, x, y, N1, N2;
	float		**f;

	read2Dfloatimage(&f, &N1, &N2, &fits, stdin);
	write2Dfloatimage(f, fits);
exit;
	fitso = copyfitsheader(fits);
	removenamedcomments("HISTORY", fitso);
	removenamedcomments("FOO_0", fitso);
	fitso->extbitpix = FLOAT_PIXTYPE;
	fitso->linebuffer = NULL;
	fitso->bscaling = 0;
	writefitsheader(fitso);
	for (y = 0; y < fits->n[1]; y++) {
		writefitsline((void *) (f[y]), fitso);	
	}
	writefitstail(fitso);
	exit(0);
}
