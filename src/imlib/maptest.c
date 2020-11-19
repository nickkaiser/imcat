/*
 * smtest.c
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "../utils/arrays.h"
#include "map.h"
#include "fits.h"

void deflect(float ri, float rj, float *di, float *dj);
int	N = 5;
int	N1, N2, M1, M2;

main(int argc, char *argv[])
{
	int	i, j, mapmode;
	float	**fsource, **ftarget;
	int	comc;
	char	*comv[MAX_COMMENTS];

	if (argc != 2) {
		fprintf(stderr, "maptest: usage: maptest mapmode\n");
		exit(-1);
	}
	sscanf(argv[1], "%d", &mapmode);
	if ((mapmode != FORWARDMAPMODE) && (mapmode != INVERSEMAPMODE)) {
		fprintf(stderr, "maptest: illegal mapmode\n");
		exit(-1);
	}
	
	fread_fits(&fsource, &M1, &M2, &comc, comv);
	N1 = M1;
	N2 = M2;
	allocFloatArray(&ftarget, N1, N2);
	if (mapmode == FORWARDMAPMODE) {
		map(ftarget, N1, N2, fsource, M1, M2, deflect);
	} else {
		set_triangle_map_mode(INVERSEMAPMODE);
		map(fsource, N1, N2, ftarget, M1, M2, deflect);
	}
	set_output_pixtype(FLOAT_PIXTYPE);	
	fwrite_fits(ftarget, N1, N2, comc, comv);	
	exit(0);
}



void deflect(float ri, float rj, float *di, float *dj)
{
/*	scrunch by 8
	*di = -0.875 * ri;
	*dj = -0.875 * rj;
*/
/*	displace
	*di = 0.1;
	*dj = 0.1;
*/
/*	stretch x 2 in y */
	*di = -0.5 * (ri - 0.5 * N2);
	*dj = 0.0;
}

