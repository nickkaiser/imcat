#define usage "\n\n\n\
NAME\n\
	flip --- reflect or rotate a fits image\n\
SYNOPSIS\n\
	flip [h|v]\n\
\n\
DESCRIPTION\n\
	By default \"flip\" rotates an image by 180 degrees to\n\
	produce an image of the same dimensions.\n\
	With optional first argument 'h' or 'v' it will reflect the image\n\
	about the horizontal or vertical axis respectively.\n\
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

#define		ROTATE_MODE 	0
#define		REFLECT_H_MODE	1
#define		REFLECT_V_MODE	2

int		main(int argc, char *argv[])	
{
	int		i, j, mode, N1, N2;
	fitsheader	*fits;
	float		**fin, **fout;
	
	mode = ROTATE_MODE;
	if (argc == 2) {
		switch (argv[1][0]) {
			case 'h':
				mode = REFLECT_H_MODE;
				break;
			case 'v':
				mode = REFLECT_V_MODE;
				break;
			default:
				error_exit(usage);
		}
	}

	read2Dfloatimage(&fin, &N1, &N2, &fits, stdin);
	allocFloatArray(&fout, N1, N2);
	
	switch (mode) {
		case ROTATE_MODE:
			for (i = 0; i < N2; i++) {
				for (j = 0; j < N1; j++) {
					fout[N2-i-1][N1-j-1] = fin[i][j];
				}
			}
			break;
		case REFLECT_H_MODE:
			for (i = 0; i < N2; i++) {
				for (j = 0; j < N1; j++) {
					fout[N2-i-1][j] = fin[i][j];
				}
			}
			break;
		case REFLECT_V_MODE:
			for (i = 0; i < N2; i++) {
				for (j = 0; j < N1; j++) {
					fout[i][N1-j-1] = fin[i][j];
				}
			}
			break;
		default:
			error_exit("flip: bad mode\n");		
	}
	
	add_comment(argc, argv, fits);
	write2Dfloatimage(fout, fits);
	exit(0);
}



