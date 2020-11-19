/*
 * fitstack_read.c
 */

#include <stdio.h>
#include <math.h>
#include "utils/error.h"
#include "fitstack_read.h"


/* number of objects per input buffer */
#define BUFF_SIZE 1000
/* number of fields */
#define N_FIELDS 8
/* max number of input buffers */
#define MAX_BUFFS 10000

/* global arrays */
double 	***x, **mag, ***r, ***z;
int	**e;

int	readmergedcat(int nexp)
{
	FILE	*lcpipe;
	int 	i, j, ibuff, nbuff, nprs, ipr, maxprs, nprsinbuff;
	double	*buff[MAX_BUFFS];

	/* open pipe */
	if (!(lcpipe = popen("lc -b -o x exp mag", "r"))) {
		fprintf(stderr, "fitstack: readmergedcat: unable to open lc-pipe for input\n");
		exit(-1);
	}

	/* read catalogue into array of buffers */
	nbuff = 0;
	nprs = nprsinbuff = 0;
	if (!(buff[nbuff] = (double *) calloc(N_FIELDS * BUFF_SIZE, sizeof(double)))) {
		error_exit("fitstack: readmergedcat: failed to allocate input buffer\n");
	}
	while (fread(buff[nbuff] + nprsinbuff * N_FIELDS, sizeof(double), 8, lcpipe)) {
		nprsinbuff++;
		nprs++;
		if (nprsinbuff == BUFF_SIZE) {	/* need a new buffer */
			if (nbuff++ >= MAX_BUFFS) {
				error_exit("fitstack: readmergedcat: too many objects in input cat\n");
			}
			if (!(buff[nbuff] = (double *) calloc(N_FIELDS * BUFF_SIZE, sizeof(double)))) {
				error_exit("fitstack: readmergedcat: failed to allocate input buffer\n");
			}
			nprsinbuff = 0;
		}
	}
	nbuff++;

	/* create data arrays */
	x = (double ***) calloc(2, sizeof(double **));
	e = (int **) calloc(2, sizeof(int *));
	mag = (double **) calloc(2, sizeof(double *));
	for (i = 0; i < 2; i++) {
		x[i] = (double **) calloc(2, sizeof(double *));
		e[i] = (int *) calloc(nprs, sizeof(int));
		mag[i] = (double *) calloc(nprs, sizeof(double));
		for (j = 0; j < 2; j++) {
			x[i][j] = (double *) calloc(nprs, sizeof(double));
		}
	}

	/* copy data from input buffers */
	nprs = 0;
	for (ibuff = 0; ibuff < nbuff; ibuff++) {
		if (ibuff == (nbuff - 1)) {
			maxprs = nprsinbuff;
		} else {
			maxprs = BUFF_SIZE;
		}
		for (ipr = 0; ipr < maxprs; ipr++) {
			x[0][0][nprs] = buff[ibuff][ipr * N_FIELDS];
			x[0][1][nprs] = buff[ibuff][ipr * N_FIELDS + 1];
			x[1][0][nprs] = buff[ibuff][ipr * N_FIELDS + 2];
			x[1][1][nprs] = buff[ibuff][ipr * N_FIELDS + 3];
			e[0][nprs] = (int) buff[ibuff][ipr * N_FIELDS + 4];
			e[1][nprs] = (int) buff[ibuff][ipr * N_FIELDS + 5];
			if (e[0][nprs] < 0 || e[0][nprs] >= nexp || e[1][nprs] < 0 || e[1][nprs] >= nexp) {
				error_exit("fitstack: readmergedcat: exposure number out of bounds\n");
			}
			mag[0][nprs] = buff[ibuff][ipr * N_FIELDS + 6];
			mag[1][nprs] = buff[ibuff][ipr * N_FIELDS + 7];
			nprs++;
		}
		free(buff[ibuff]);
	}
	return (nprs);
}  
