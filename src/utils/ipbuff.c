/*
 * ipbuff.c
 */


#define	BUFFPTRSIZE	10000
#define	BUFFSIZE	10000

#include <stdio.h>
#include <stdlib.h>
#include "ipbuff.h"

double	**readdoublebuff(int width, FILE *stream, int *np)
{
	double	**buffp,		/* pointer returned */
		*buffptr[BUFFPTRSIZE],	/* array of pointers to input buffers */
		*buff;			/* input buffer */
	int	ibuff, p, pp;		/* indices */

	/* allocate 1st input buffer array */
	ibuff = 0;
	buffptr[ibuff] = buff = (double *) calloc(width * BUFFSIZE, sizeof(double));
	p = pp = 0;
	while (fread(buff + pp * width, sizeof(double), width, stream)) {
		p++;
		pp++;
		if (pp == BUFFSIZE) {
			ibuff++;
			if (ibuff >= BUFFPTRSIZE) {
				error_exit("catipbuff: too many objects in input cat\n");
			}
			buffptr[ibuff] = buff = (double *) calloc(width * BUFFSIZE, sizeof(double));
			pp = 0;
		}
	}
	*np = p;		/* set total number of particles */

	/* allocate returned array and assign pointers in it */
	buffp = (double **) calloc(*np, sizeof(double));
	while (ibuff >= 0) {
		p--;
		pp--;
		if (pp < 0) {
			ibuff--;
			if (ibuff < 0) {
				break;
			}
			pp = BUFFSIZE - 1;
		}
		buffp[p] = buffptr[ibuff] + pp * width;
	}

	return(buffp);
}