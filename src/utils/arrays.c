/**
 ** routines for allocating, freeing and copying 2-D/3-D arrays
 **
 ** these are the newstyle contiguous storage versions
 **
 **/
 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "error.h"
#include "arrays.h"

void	*checkalloc(int nel, int size)
{
	void *data;
	data = calloc(nel, size);
	if (!data) {
		error_exit("checkalloc: memory allocation failure\n");
	}
	return (data);
}


void	allocShortArray(short ***f, int N1, int N2)
{
	int		i;
	
	(*f) = (short **) calloc(N2, sizeof(short *));
	if (!*f)
		error_exit("allocShortArray: memory allocation failure\n");
	(*f)[0] = (short *) calloc(N1 * N2, sizeof(short));
	if (!(*f)[0])
		error_exit("allocShortArray: memory allocation failure\n");
	for (i = 1; i < N2; i++)
		(*f)[i] = (*f)[i - 1] + N1;
}




void	allocFloatArray(float ***f, int N1, int N2)
{
	int		i;
	
	(*f) = (float **) calloc(N2, sizeof(float *));
	if (!*f)
		error_exit("allocFloatArray: memory allocation failure\n");
	(*f)[0] = (float *) calloc(N1 * N2, sizeof(float));
	if (!(*f)[0])
		error_exit("allocFloatArray: memory allocation failure\n");
	for (i = 1; i < N2; i++)
		(*f)[i] = (*f)[i - 1] + N1;
}





void	freeShortArray(short **f, int N1, int N2)
{
	int	i;
	
	free(f[0]);
	free(f);
}





void	freeFloatArray(float **f, int N1, int N2)
{
	int	i;
	
	free(f[0]);
	free(f);
}





void	copyFloatToShort(float **fsrc, short **fdst, int N1, int N2)
{
	int	i, j;
	
	for (i = 0; i < N2; i++)
		for (j = 0; j < N1; j++)
			fdst[i][j] = fsrc[i][j];
}





void	copyShortToFloat(short **fsrc, float **fdst, int N1, int N2)
{
	int	i, j;
	
	for (i = 0; i < N2; i++)
		for (j = 0; j < N1; j++)
			fdst[i][j] = fsrc[i][j];
}



float	***alloc3DFloatArray(int N1, int N2, int N3)
{
	int		i, j, nel;
	float		***f, *fdat;
	
	nel = N1 * N2 * N3;
	fdat = (float *) checkalloc(nel, sizeof(float));
	f = (float ***) checkalloc(N3, sizeof(float));
	for (i = 0; i < N3; i++) {
		f[i] = (float **) checkalloc(N2, sizeof(float));
		for (j = 0; j < N2; j++) {
			f[i][j] = fdat;
			fdat += N1;
		}
	}
	return (f);
}

void	free3DFloatArray(float ***f, int N1, int N2, int N3)
{
	int	i;

	free(f[0][0]);
	for (i = 0; i < N3; i++) {
		free(f[i]);
	}
	free(f);
}
