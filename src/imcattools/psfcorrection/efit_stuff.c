/*
 * efit_stuff.c
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../../utils/error.h"
#include "efit_stuff.h"

static int	framesize;

double		g(int mode, double xx, double yy)
{
	double	x[2];
	
	x[0] = (xx - framesize / 2);
	x[1] = (yy - framesize / 2);
	
	switch (mode) {
		case 0:
			return(1.0);
			break;
		case 1:
			return(x[0]);
			break;
		case 2:
			return(x[1]);
			break;
		case 3:
			return(x[0] * x[0]);
			break;
		case 4:
			return(x[0] * x[1]);
			break;
		case 5:
			return(x[1] * x[1]);
			break;
		case 6:
			return(x[0] * x[0] * x[0]);
			break;
		case 7:
			return(x[0] * x[0] * x[1]);
			break;
		case 8:
			return(x[0] * x[1] * x[1]);
			break;
		case 9:
			return(x[1] * x[1] * x[1]);
			break;
		case 10:
			return(x[0] * x[0] * x[0] * x[0]);
			break;
		case 11:
			return(x[0] * x[0] * x[0] * x[1]);
			break;
		case 12:
			return(x[0] * x[0] * x[1] * x[1]);
			break;
		case 13:
			return(x[0] * x[1] * x[1] * x[1]);
			break;
		case 14:
			return(x[1] * x[1] * x[1] * x[1]);
			break;
		case 15:
			return(x[0] * x[0] * x[0] * x[0] * x[0]);
			break;
		case 16:
			return(x[0] * x[0] * x[0] * x[0] * x[1]);
			break;
		case 17:
			return(x[0] * x[0] * x[0] * x[1] * x[1]);
			break;
		case 18:
			return(x[0] * x[0] * x[1] * x[1] * x[1]);
			break;
		case 19:
			return(x[0] * x[1] * x[1] * x[1] * x[1]);
			break;
		case 20:
			return(x[1] * x[1] * x[1] * x[1] * x[1]);
			break;
		case 21:
			return(x[0] * x[0] * x[0] * x[0] * x[0] * x[0]);
			break;
		case 22:
			return(x[0] * x[0] * x[0] * x[0] * x[0] * x[1]);
			break;
		case 23:
			return(x[0] * x[0] * x[0] * x[0] * x[1] * x[1]);
			break;
		case 24:
			return(x[0] * x[0] * x[0] * x[1] * x[1] * x[1]);
			break;
		case 25:
			return(x[0] * x[0] * x[1] * x[1] * x[1] * x[1]);
			break;
		case 26:
			return(x[0] * x[1] * x[1] * x[1] * x[1] * x[1]);
			break;
		case 27:
			return(x[1] * x[1] * x[1] * x[1] * x[1] * x[1]);
			break;
		default:
			error_exit("g: bad mode\n");
			break;
	}
}





void		readamplitudes(int *nmodes, int *framesize, int *order, double ***amp, FILE *stream)
{
	char	line[1024];
	int		mode, junk;
	
	if (!fgets(line, 1024, stream))
		error_exit("readamplitudes: can't read the amplitudes file\n");
	if (!sscanf(line, "%d", nmodes))
		error_exit("readamplitudes: can't read nmodes from amplitudes file\n");
	if (!fgets(line, 1024, stream))
		error_exit("readamplitudes: can't read the amplitudes file\n");
	if (!sscanf(line, "%d", order))
		error_exit("readamplitudes: can't read order from amplitudes file\n");
	if (!fgets(line, 1024, stream))
		error_exit("readamplitudes: can't read the amplitudes file\n");
	if (!sscanf(line, "%d", framesize))
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
		if (3 != sscanf(line, "%d %lf %lf", &junk, (*amp)[0]+mode, (*amp)[1]+mode))
			error_exit("readamplitudes: can't read amplitude from amplitudes file\n");
		mode++;
	}
	if (junk != *nmodes - 1)
		error_exit("readamplitudes: number of last mode read doesn't match nmodes - 1\n");
}


void		writeamplitudes(int nmodes, int order, int framesize, double **amp)
{
		int	mode;

		fprintf(stdout, "%d modes\n", nmodes);
		fprintf(stdout, "%d order model\n", order);
		fprintf(stdout, "%d pixels wide\n", framesize);
		for (mode = 0; mode < nmodes; mode++) {
			fprintf(stdout, "%d %13.8lg %13.8lg\n", mode, amp[0][mode], amp[1][mode]);
		}
}

void		setframesize(int N)
{
	framesize = N;
}



