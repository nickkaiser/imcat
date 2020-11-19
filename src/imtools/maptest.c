/*
 * smtest.c
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "../utils/arrays.h"
#include "../imlib/map.h"

void deflect(float ri, float rj, float *di, float *dj);
int	N = 5;

main(int argc, char *argv[])
{
	int		i, j;
	float	**fsource, **ftarget;
	
	allocFloatArray(&fsource, N, N);
	allocFloatArray(&ftarget, N, N);
	
	fprintf(stderr, "source image: \n");
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			fsource[i][j] = N * i + j;
			fprintf(stderr, "%10.3f ", fsource[i][j]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
	
	map(ftarget, N, N, fsource, N, N, deflect);
	
	fprintf(stderr, "target image: \n");
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			fsource[i][j] = N * i + j;
			fprintf(stderr, "%10.3f ", ftarget[i][j]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "used %ld triangles in total\n", trianglecount());
	exit(0);
}



void deflect(float ri, float rj, float *di, float *dj)
{
	*di = N - 2 * ri + 0.1;
	*dj = N - 2 * rj - 0.1;
}

