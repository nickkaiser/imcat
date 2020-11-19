/*
 * fits.c
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../utils/error.h"
#include "fits.h"
#include "../utils/arrays.h"

static	FILE	*theipf = stdin, *theopf = stdout;

/* default pixel types */
static	int	bitpixout = SHORT_PIXTYPE;
static	int	bitpixin  = SHORT_PIXTYPE;

void	set_output_pixtype(int thetype)
{
	switch (thetype) {
		case SHORT_PIXTYPE:
			bitpixout = SHORT_PIXTYPE;
			break;
		case FLOAT_PIXTYPE:
			bitpixout = FLOAT_PIXTYPE;
			break;
		case INT_PIXTYPE:
			bitpixout = INT_PIXTYPE;
			break;
		case DBL_PIXTYPE:
			bitpixout = FLOAT_PIXTYPE;
			break;
		default:
			error_exit("set_output_pixtype: bad pixtype\n");
			break;
	}
}


void	set_input_pixtype(int thetype)
{
	switch (thetype) {
		case SHORT_PIXTYPE:
			bitpixin = SHORT_PIXTYPE;
			break;
		case FLOAT_PIXTYPE:
			bitpixin = FLOAT_PIXTYPE;
			break;
		case INT_PIXTYPE:
			bitpixin = INT_PIXTYPE;
			break;
		case DBL_PIXTYPE:
			bitpixin = DBL_PIXTYPE;
			break;
		default:
			error_exit("set_output_pixtype: bad pixtype\n");
			break;
	}
}

void	get_input_pixtype(int *thetype)
{
	*thetype = bitpixin;
}



void	write_fits_head(int N1, int N2, int comc, char *comv[])
{
	int 	com = 0, icom = 0, ic;
	char	record[2880], line[COM_LENGTH];
	long	length;

	for (ic = 0; ic < 2880; ic++)
		record[ic] = ' ';
	sprintf(line, "SIMPLE  =                    T / written by imcat");
	strncpy(record + icom++, line, 49);
	sprintf(line, "BITPIX  =  %19d /", bitpixout);
	strncpy(record + icom++ * COM_LENGTH, line, 32);
	sprintf(line, "NAXIS   =                    2 /");
	strncpy(record + icom++ * COM_LENGTH, line, 32);
	sprintf(line, "NAXIS1  =  %19d /", N1);
	strncpy(record + icom++ * COM_LENGTH, line, 32);
	sprintf(line, "NAXIS2  =  %19d /", N2);
	strncpy(record + icom++ * COM_LENGTH, line, 32);
	for (com = 5; com < comc; com++) {
		strncpy(record + icom++ * COM_LENGTH, comv[com], COM_LENGTH);
		if (icom == 36) {
			fwrite(record, sizeof(char), 2880, theopf);	
			for (ic = 0; ic < 2880; ic++)
				record[ic] = ' ';
			icom = 0;
		}
	}	
	sprintf(line, "END");
	strncpy(record + icom++ * COM_LENGTH, line, 3);
	fwrite(record, sizeof(char), 2880, theopf);	
}


void	write_fits_line(short *f, int N1)
{
	int i;

	fwrite(f, sizeof(short), N1, theopf);
}

void	fwrite_fits_line(float *f, int N1)
{
	int i;
	short	*fshort;
	int	*fint, fi;
	double	*fdouble;

	switch (bitpixout) {
		case SHORT_PIXTYPE:
			fshort = (short *) f;
			for (i = 0; i < N1; i++) {
				fi = (int) floor(0.5 + f[i]);
				fshort[i] = (short) (fi > SHRT_MAX ? SHRT_MAX : (fi < SHRT_MIN ? (SHRT_MIN + 1): fi));
			}
			fwrite(fshort, sizeof(short), N1, theopf);
			break;
		case FLOAT_PIXTYPE:
			fwrite(f, sizeof(float), N1, theopf);
			break;
		case INT_PIXTYPE:
			fint = (int *) f;
			for (i = 0; i < N1; i++)
				fint[i] = (int) floor(0.5 + f[i]);
			fwrite(fint, sizeof(int), N1, theopf);
			break;
		default:
			error_exit("fwrite_fits_line: bad pixtype\n");
	}
}

void	write_fits_tail(int N1, int N2)
{
	int 	ic, pad;
	char	record[2880];
	long	length;

	switch (bitpixout) {
		case SHORT_PIXTYPE:
			length = (long) N1 * (long) N2 * sizeof(short);
			break;
		case FLOAT_PIXTYPE:
			length = (long) N1 * (long) N2 * sizeof(float);
			break;
		case INT_PIXTYPE:
			length = (long) N1 * (long) N2 * sizeof(int);
			break;
		default:
			error_exit("write_fits_tail: illegal pixtype\n");
			break;
	}
	pad = 2880 - length + 2880 * (length / 2880);
	for (ic = 0; ic < pad; ic++)
		record[ic] = ' ';
	fwrite(record, sizeof(char), pad, theopf);	
}


void	write_fits(short **f, int N1, int N2, int comc, char *comv[])
{
	int	i;

	write_fits_head(N1, N2, comc, comv);
	for (i = 0; i < N2; i++)
		write_fits_line(f[i], N1);
	write_fits_tail(N1, N2);
}



void	fwrite_fits(float **f, int N1, int N2, int comc, char *comv[])
{
	int	i, ix, iy;
	short	*fshort;
	float	*ffloat, fi;
	int	*fint;

	write_fits_head(N1, N2, comc, comv);

	switch (bitpixout) {
		case SHORT_PIXTYPE:
			fshort = (short *) calloc(N1 * N2, sizeof(short));
			i = 0;
			for (iy = 0; iy < N2; iy++) {
				for (ix = 0; ix < N1; ix++) {
					fi = f[iy][ix];
					fshort[i++] = (short) (fi > SHRT_MAX ? SHRT_MAX : (fi < SHRT_MIN ? (SHRT_MIN + 1): fi));
				}
			}
			fwrite(fshort, sizeof(short), N1 * N2, theopf);
			free(fshort);
			break;
		case FLOAT_PIXTYPE:
			ffloat = (float *) calloc(N1 * N2, sizeof(float));
			i = 0;
			for (iy = 0; iy < N2; iy++) {
				for (ix = 0; ix < N1; ix++) {
					ffloat[i++] = f[iy][ix];
				}
			}
			fwrite(ffloat, sizeof(float), N1 * N2, theopf);
			free(ffloat);
			break;
		case INT_PIXTYPE:
			fint = (int *) calloc(N1 * N2, sizeof(int));
			i = 0;
			for (iy = 0; iy < N2; iy++) {
				for (ix = 0; ix < N1; ix++) {
					fint[i++] = (int) floor(0.5 + f[iy][ix]);
				}
			}
			fwrite(fint, sizeof(int), N1 * N2, theopf);
			free(fint);
			break;
		default:
			error_exit("fwrite_fits: bad pixtype\n");	
	}
	write_fits_tail(N1, N2);
}



void	read_fits_head(int *N1, int *N2, int *comc, char *comv[])
{
	int 	i, naxis, bitpix, lineno, done = 0;
	char	record[2880];	
	char	*line;
	
	fread(record, sizeof(char), 2880, theipf);
	*comc = 0;
	for (lineno = 0; lineno < 5; lineno++) {
		line = record + lineno * COM_LENGTH;
		comv[*comc] = (char *) calloc(COM_LENGTH, sizeof(char));
		strncpy(comv[(*comc)++], line, COM_LENGTH);
		switch (lineno) {
			case 0:
				break;
			case 1:
				sscanf(line + 9, "%d", &bitpix);
				switch (bitpix) {
					case SHORT_PIXTYPE:
					case FLOAT_PIXTYPE:
					case INT_PIXTYPE:
					case DBL_PIXTYPE:
						bitpixin = bitpix;
						break;
					default:
						error_exit("read_fits_head: illegal pixel type\n");
						break; 
				}
				break;
			case 2:
				sscanf(line + 9, "%d", &naxis);
				break;
			case 3:
				sscanf(line + 9, "%d", N1);
				break;
			case 4:
				sscanf(line + 9, "%d", N2);
				break;
		}
	}	
	lineno = 5;
	while (!done) {
		for (; lineno < 36 ; lineno++) {
			line = record + lineno * COM_LENGTH;
			if (!strncmp(line, "END", 3)) {
				done = 1;
				break;
			} 
/*			if (    !strncmp(line, "EXPNUM  =", 9) ||
				!strncmp(line, "OBJECT  =", 9) ||
				!strncmp(line, "RA      =", 9) ||
				!strncmp(line, "DEC     =", 9) ||
				!strncmp(line, "DATEOBS =", 9) ||
				!strncmp(line, "INTTIMER=", 9) ||
				!strncmp(line, "INTTIME =", 9) ||
				!strncmp(line, "UTIME   =", 9) ||
				!strncmp(line, "DETECTOR=", 9) ||
				!strncmp(line, "INSTRUME=", 9) ||
				!strncmp(line, "BSCALE  =", 9) ||
				!strncmp(line, "BZERO   =", 9) ||
				!strncmp(line, "ESCALE  =", 9) ||
				!strncmp(line, "PIXSIZE =", 9) ||
				!strncmp(line, "AIRMASS =", 9) ||
				!strncmp(line, "FILTER  =", 9) ||
				!strncmp(line, "OBSTYPE =", 9) ||
				!strncmp(line, "HISTORY", 7)) */
			{
				comv[*comc] = (char *) calloc(COM_LENGTH, sizeof(char));
				strncpy(comv[*comc], line, COM_LENGTH);
				(*comc)++;
			} 
		}
		if (!done) 
			fread(record, sizeof(char), 2880, theipf);
		lineno = 0;	
	}
}

void	read_fits_line(short *f, int N1)
{
	int	i;

	fread(f, sizeof(short), N1, theipf);
}

void	fread_fits_line(float *f, int N1)
{
	int	i;
	short	*fshort;
	int	*fint;
	double	*fdouble;

	switch (bitpixin) {
		case DBL_PIXTYPE:
			/* totally inefficient kluge - hopefully we rarely use it */
			fdouble = (double *) calloc(N1, sizeof(double));
			fread(fdouble, sizeof(double), N1, theipf);
			for (i = N1 - 1; i >= 0; i--)
				f[i] = (float) fdouble[i];
			free(fdouble);
			break;
		case SHORT_PIXTYPE:
			fshort = (short *) f;
			read_fits_line(fshort, N1);
			for (i = N1 - 1; i >= 0; i--)
				f[i] = (float) fshort[i];
			break;
		case FLOAT_PIXTYPE:
			fread(f, sizeof(float), N1, theipf);
			fint = (int *) f;
			break;
		case INT_PIXTYPE:
			fint = (int *) f;
			fread(fint, sizeof(int), N1, theipf);
			for (i = 0; i < N1; i++)
				f[i] = (float) fint[i];
			break;
		default:
			error_exit("fwrite_fits_line: bad pixtype\n");
	}

}

void	read_fits(short ***f, int *N1, int *N2, int *comc, char *comv[])
{
	int 	i;
	
	read_fits_head(N1, N2, comc, comv);
	
	*f = (short **) calloc(*N2, sizeof(short *));
	if (!*f) error_exit("read_fits: memory allocation error");
	for (i = 0; i < *N2; i++)
		{
			(*f)[i] = (short *) calloc(*N1, sizeof(short));
			if (!(*f)[i]) 
				error_exit("read_fits: memory allocation error");
		}
	for (i = 0; i < *N2; i++) {
		read_fits_line((*f)[i], *N1);
	}
}




void	fread_fits(float ***f, int *N1, int *N2, int *comc, char *comv[])
{
	int 	i, j;
	
	read_fits_head(N1, N2, comc, comv);
	allocFloatArray(f, *N1, *N2);
	for (i = 0; i < *N2; i++) {
		fread_fits_line((*f)[i], *N1);
	}
}






void	argsToString(int argc, char **argv, char *string)
{
	int	i, j, pos = 0, len, left = COM_LENGTH;
	
	for (i = 0; i < argc; i++) {
		len = strlen(argv[i]);
		for (j = 0; j < len; j++) {
			string[pos++] = argv[i][j];
			if (pos == COM_LENGTH - 2)
				break;
		}
		string[pos++] = ' ';
		if (pos > COM_LENGTH - 2)
			break;
	}
	string[pos] = '\0';
}




void	add_comment(int argc, char **argv, int *comc, char **comv)
{
	char	argstring[COM_LENGTH];

	if (*comc > MAX_COMMENTS - 1)
		return;
	if (*comc < 5)
		*comc = 5;

	comv[*comc] = (char *) calloc(COM_LENGTH, sizeof(char));
	argsToString(argc, argv, argstring);
	sprintf(comv[*comc], "HISTORY = ");
	strncpy(comv[*comc] + 10, argstring, COM_LENGTH - 12);
	comv[*comc][COM_LENGTH-1] = '\0';
	(*comc)++;
}


void	add_1comment(char *comment, int *comc, char **comv)
{
	if (*comc > MAX_COMMENTS - 1)
		return;

	comv[*comc] = (char *) calloc(COM_LENGTH, sizeof(char));
	sprintf(comv[*comc], "HISTORY = ");
	strncpy(comv[*comc] + 10, comment, COM_LENGTH - 12);
	comv[*comc][COM_LENGTH-1] = '\0';
	(*comc)++;
}



void	addheadervalue(char *name, char *value, int *comc, char **comv)
{
	if (*comc > MAX_COMMENTS - 1)
		error_exit("addheadervalue: too many header lines\n");
	if (strlen(name) > 8)
		error_exit("addheadervalue: name too long\n");
	if (strlen(name) > 68)
		error_exit("addheadervalue: value too long\n");
	comv[*comc] = (char *) calloc(COM_LENGTH, sizeof(char));
	sprintf(comv[*comc], "%-8s= %-68s", name, value);
	comv[*comc][COM_LENGTH-1] = '\0';
	(*comc)++;
}



char	*getheadervalue(char *name, int comc, char **comv)
{
	int	com, thecom = -1;
	char 	*value;

	for (com = 0; com < comc; com++) {
		if (!strncmp(name, comv[com], 8)) {
			thecom = com;
		}
	}
	if (thecom < 0)
		return (NULL);
	value = (char *) calloc(1 + strlen(comv[thecom] + 10), sizeof(char));
	strcpy(value, comv[thecom] + 10);
	return (value);
}



void	set_fits_opf(FILE *opf)
{
	theopf = opf;
}



void	set_fits_ipf(FILE *ipf)
{
	theipf = ipf;
}


/*
 * multidimensional image routines
 */

void	read_fits_head_ND(int *Naxes, int *N, int *comc, char *comv[])
{
	int 	i, naxes, bitpix, lineno, done = 0;
	char	record[2880];	
	char	*line;
	
	fread(record, sizeof(char), 2880, theipf);
	sscanf(record + 1 * COM_LENGTH + 9, "%d", &bitpix);
	switch (bitpix) {
		case SHORT_PIXTYPE:
		case FLOAT_PIXTYPE:
		case INT_PIXTYPE:
		case DBL_PIXTYPE:
			bitpixin = bitpix;
			break;
		default:
			error_exit("read_fits_head: illegal pixel type\n");
			break; 
	}
	sscanf(record + 2 * COM_LENGTH + 9, "%d", &naxes);
	*Naxes = naxes;
	if (naxes > 7)
		error_exit("read_fits_head_ND: too many axes\n");
	for (i = 0; i < naxes; i++)
		sscanf(record + (3 + i) * COM_LENGTH + 9, "%d", &(N[i]));
	*comc = 0;
	lineno = 3 + naxes;
	while (!done) {
		for (; lineno < 36 ; lineno++) {
			line = record + lineno * COM_LENGTH;
			if (!strncmp(line, "END", 3)) {
				done = 1;
				break;
			} 
			if (!strncmp(line, "OBJECT  =", 9) || 
				!strncmp(line, "HISTORY =", 9)) {
				comv[*comc] = (char *) calloc(COM_LENGTH, sizeof(char));
				strncpy(comv[*comc], line, COM_LENGTH);
				(*comc)++;
			} 
		}
		if (!done) 
			fread(record, sizeof(char), 2880, theipf);
		lineno = 0;	
	}
}




void	write_fits_head_ND(int Naxes, int *N, int comc, char *comv[])
{
	int 	i, com = 0, icom, ic;
	char	record[2880], line[COM_LENGTH];
	long	length;
	
	for (ic = 0; ic < 2880; ic++)
		record[ic] = ' ';
	sprintf(line, "SIMPLE  =                    T / written by imcat");
	strncpy(record, line, 49);
	sprintf(line, "BITPIX  =  %19d /", bitpixout);
	strncpy(record + 1 * COM_LENGTH, line, 32);
	sprintf(line, "NAXIS   =  %19d /", Naxes);
	strncpy(record + 2 * COM_LENGTH, line, 32);
	for (i = 0; i < Naxes; i++) {
		sprintf(line, "NAXIS%1d  =  %19d /", i + 1, N[i]);
		strncpy(record + (3 + i) * COM_LENGTH, line, 32);
	}
	icom = 3 + Naxes;
	for (com = icom; com < comc; com++) {
		strncpy(record + icom++ * COM_LENGTH, comv[com], COM_LENGTH);
		if (icom == 36) {
			fwrite(record, sizeof(char), 2880, theopf);	
			for (ic = 0; ic < 2880; ic++)
				record[ic] = ' ';
			icom = 0;
		}
	}
	sprintf(line, "END");
	strncpy(record + icom++ * COM_LENGTH, line, 3);
	fwrite(record, sizeof(char), 2880, theopf);	
}



void	write_fits_tail_ND(int Naxes, int *N)
{
	int 	i, ic, pad;
	char	record[2880];
	long	length;

	switch (bitpixout) {
		case SHORT_PIXTYPE:
			length = (long) sizeof(short);
			break;
		case FLOAT_PIXTYPE:
			length = (long) sizeof(float);
			break;
		case INT_PIXTYPE:
			length = (long) sizeof(int);
			break;
		default:
			error_exit("write_fits_tail: illegal pixtype\n");
			break;
	}
	for (i = 0; i < Naxes; i++)
		length *= (long) N[i];
	pad = 2880 - length + 2880 * (length / 2880);
	for (ic = 0; ic < pad; ic++)
		record[ic] = ' ';
	fwrite(record, sizeof(char), pad, theopf);	
}




