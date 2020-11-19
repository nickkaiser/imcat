#define usage "\n\n\n\
SYNOPSIS\n\
	addfitscomment fitsfile name value\n\
\n\
DESCRIPTION\n\
	\"addfitscomment\" appends a comment to the fits HISTORY list.\n\
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
#include "../utils/arrays.h"

int		main(int argc, char *argv[])	
{
/*
	int		N1, N2;
	int		comc, pixtype;
	char	*comv[MAX_COMMENTS];
	float	**f;
	char	filename, comment;
	FILE	*fitsf;
*/
	error_exit("addfitscomment: has been superseded by ic\n");

	/*
	if (argc != 4)
		error_exit(usage);
	if (strlen(argv[2]) > 8)
		error_exit("addfitscomment: name too long\n");
	if (strlen(argv[3]) > 68)
		error_exit("addfitscomment: name too long\n");
	
	fitsf = fopen(argv[1], "r");
	if (!fitsf)
		error_exit("addfitscomment: failed to open fits file for input\n");
	set_fits_ipf(fitsf);
	fread_fits(&f, &N1, &N2, &comc, comv);
	get_input_pixtype(&pixtype);
	fclose(fitsf);
	addheadervalue(argv[2], argv[3], &comc, comv);
	fitsf = fopen(argv[1], "w");
	if (!fitsf)
		error_exit("addfitscomment: failed to open fits file for output\n");
	set_fits_opf(fitsf);
	set_output_pixtype(pixtype);
	fwrite_fits(f, N1, N2, comc, comv);
	fclose(fitsf);
	exit(0);
	*/
}

