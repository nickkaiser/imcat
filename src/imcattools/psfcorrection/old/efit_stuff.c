/*
 * efit_stuff.c
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../../utils/error.h"
#include "efit_stuff.h"

static int	nframes, framesize;

double		g(int mode, int frame, double xx, double yy)
{
	double	x[2];
	
	x[0] = (xx - framesize / 2);
	x[1] = (yy - framesize / 2);
	
	if (mode < nframes)
		return ((mode == frame ? 1.0 : 0.0));
	switch (mode - nframes) {
		case GRAD1:
			return(x[0]);
			break;
		case GRAD2:
			return(x[1]);
			break;
		case GRADGRAD11:
			return(x[0] * x[0]);
			break;
		case GRADGRAD12:
			return(x[0] * x[1]);
			break;
		case GRADGRAD22:
			return(x[1] * x[1]);
			break;
		default:
			error_exit("g: bad mode\n");
			break;
	}
}





void		readamplitudes(int *nmodes, int *nframes, int *framesize, 
			int *order, double ***amp, FILE *stream)
{
	char	line[1024];
	int		mode, junk;
	
	if (!fgets(line, 1024, stream))
		error_exit("readamplitudes: can't read the amplitudes file\n");
	if (EOF == sscanf(line, "%d", nmodes))
		error_exit("readamplitudes: can't read nmodes from amplitudes file\n");
	if (!fgets(line, 1024, stream))
		error_exit("readamplitudes: can't read the amplitudes file\n");
	if (EOF == sscanf(line, "%d", order))
		error_exit("readamplitudes: can't read order from amplitudes file\n");
	if (!fgets(line, 1024, stream))
		error_exit("readamplitudes: can't read the amplitudes file\n");
	if (EOF == sscanf(line, "%d", nframes))
		error_exit("readamplitudes: can't read nframes from amplitudes file\n");
	if (!fgets(line, 1024, stream))
		error_exit("readamplitudes: can't read the amplitudes file\n");
	if (EOF == sscanf(line, "%d", framesize))
		error_exit("readamplitudes: can't read framesize from amplitudes file\n");
	*amp = (double **) calloc(2, sizeof(double *));
	(*amp)[0] = (double *) calloc(*nmodes, sizeof(double));
	(*amp)[1] = (double *) calloc(*nmodes, sizeof(double));	
	if (!((*amp)[0]) || !((*amp)[1]) || !(*amp))
		error_exit("readamplitudes: memory allocation failed\n");
	mode = 0;
	while (fgets(line, 1024, stream)) {
		if (line[0] == '#')
			continue;
		if (EOF == sscanf(line, "%d %lf %lf", &junk, (*amp)[0]+mode, (*amp)[1]+mode))
			error_exit("readamplitudes: can't read amplitude from amplitudes file\n");
		mode++;
	}
	if (junk != *nmodes - 1)
		error_exit("readamplitudes: number of last mode read doesn't match nmodes - 1\n");
}


void		setnframes(int xnframes)
{
	nframes = xnframes;
}



void		setframesize(int N)
{
	framesize = N;
}



