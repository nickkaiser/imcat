/*
 * modefunc.c
 */

#include <stdio.h>
#include <math.h>
#include "error.h"
#include "modefunc.h"

/* origin of coordinates */
static 	double	x0[2];

/* args */
static	char	*argstring = NULL;

double	f(int l, int m, double *x)
{
	double	res = 1, xx[2];
	int	i;

	xx[0] = x[0] - x0[0];
	xx[1] = x[1] - x0[1];
	for (i = 0; i < (l - m) ; i++) {
		res *= xx[0];
	}
	for (i = 0; i < m ; i++) {
		res *= xx[1];
	}
	return (res);
}

int	write2Dpolymodel(char *parfile, int nmodes, int *l, int *m, int asize, double **a, int nvar, char *vardef[], char *xvar)
{
	int	M, ivar, i;
	FILE	*parf = stdout;
	char	lcstring[1024];
	
	sprintf(lcstring, "lc -C -x -H 'model_type = {2Dpoly}' -H 'dim = %d' -H 'nmodes = %d' -H 'xname = {%s}' -H 'xorigin = %lf %lf 2 vector' -H 'nvar = %d'",
		asize, nmodes, xvar, x0[0], x0[1], nvar);
	strcat(lcstring, " -n l -n m");
		for (ivar = 0; ivar < nvar; ivar++) {
		strcat(lcstring, " -N '");
		strcat(lcstring, vardef[ivar]);
		strcat(lcstring, "'");
	}
	strcat(lcstring, " -a 'spatial polynomial parameter file'");

	if (argstring) {
		strcat(lcstring, " -a 'history: ");
		strcat(lcstring, argstring);
		strcat(lcstring, "' ");
	}
		
		
	if (parfile) {
		strcat(lcstring, " > ");
		strcat(lcstring, parfile);
	}	
	
	parf = popen(lcstring, "w");
	if (!parf) {
		error_exit("write2Dpolymodel: failed to open lc-pipe for output\n");
	}
		
	for (M = 0; M < nmodes; M++) {
		fprintf(parf, "%3d %3d ", l[M], m[M]);
		for (i = 0; i < asize; i++) {
			fprintf(parf, "%13.8lg ", a[i][M]);
		}
		fprintf(parf, "\n");
	}
	pclose(parf);
}

void	get2Dpolymodel(char *filename, int **lptr, int**mptr, int *asize, double ***aptr, int *nmodesptr, int *nvar, char *vardef[], char **xname)
{
	int	lcformat;
	char	sysstring[128], defxname[2] = "x";

	if (filename) {
		sprintf(sysstring, "grep \"catalogue file\" %s > /dev/null", filename);
		lcformat = !system(sysstring);
	} else {
		lcformat = 1;
	}

	if (lcformat) {
		getmodeamplitudes_lc(filename, lptr, mptr, asize, aptr, nmodesptr, nvar, vardef, xname);
	} else {
		getmodeamplitudes_txt(filename, lptr, mptr, aptr, nmodesptr);
		*asize = 2;
		*nvar = 1;
		vardef[0] = "a 1 2";
		*xname = defxname;
	}
}

void	getmodeamplitudes_lc(char *filename, int **lptr, int**mptr, int *asize, double ***aptr, int *nmodesptr, int *nvar, char *vardef[], char **xname)
{
	FILE *lcpipe;
	char	lccommand[1024], headerline[1024], model_type[128], thexname[128];
	int	*l, *m, imode, nmodes, var, dim, i;
	double	**a;
	double	*inputrecord, *ld, *md;

	/* generate the lc-command */
	sprintf(lccommand, "lc -b -P model_type -P dim -P nmodes -P xname -P xorigin l m +all");
	if (filename) {
		strcat(lccommand, " < ");
		strcat(lccommand, filename);
	}
	lcpipe = popen(lccommand, "r");
	if (!lcpipe) {
		error_exit("getmodeamplitudes_lc: failed to open lc-pipe for input\n");
	}
	
	/* read the header vals */
	fgets(headerline, 1024, lcpipe);
	sscanf(headerline, "%s", model_type);
	if (strcmp(model_type, "2Dpoly")) {
		error_exit("getmodeamplitudes_lc: invalid model_type header value\n");
	}
	fgets(headerline, 1024, lcpipe);
	sscanf(headerline, "%d", &dim);
	fgets(headerline, 1024, lcpipe);
	sscanf(headerline, "%d", nmodesptr);
	nmodes = *nmodesptr;
	fgets(headerline, 1024, lcpipe);
	sscanf(headerline, "%s", thexname);
	*xname = thexname;
	fgets(headerline, 1024, lcpipe);
	sscanf(headerline, "%lf %lf", x0, x0+1);

	if (!getvars(lcpipe, nvar, vardef, asize)) {
		error_exit("getmodeamplitudes_lc: failed to read header info\n");
	}
	*asize -= 2;
	if (*asize != dim) {
		error_exit("getmodeamplitudes_lc: mismatch betwene header 'dim' value and size of object items\n");
	}
	if (strcmp(vardef[0], "1 1 l") || strcmp(vardef[1], "1 1 m")) {
		error_exit("getmodeamplitudes_lc: l, m items seem to have wrong dimension!\n");
	}

	l = (int *) calloc(nmodes, sizeof(int));
	m = (int *) calloc(nmodes, sizeof(int));
	*aptr = a = (double **) calloc(dim, sizeof(double *));
	for (i = 0; i < dim; i++) {
		a[i] = (double *) calloc(nmodes, sizeof(double));
	}
	inputrecord = (double *) calloc(dim + 2, sizeof(double));
	ld = inputrecord;
	md = inputrecord + 1;
	imode = 0;
	while (fread(inputrecord, sizeof(double), dim + 2, lcpipe)) {
		if (imode == nmodes) {
			error_exit("getmodeamplitudes_lc: parameter file is too long!\n");
		}
		l[imode] = (int) floor(inputrecord[0]);
		m[imode] = (int) floor(*md);
		for (i = 0; i < dim; i++) {
			a[i][imode] = (double) (inputrecord[2 + i]);
		}
		imode++;
	}

	*lptr = l;
	*mptr = m;
	*nmodesptr = nmodes;
}

void	getmodeamplitudes_txt(char *filename, int **lptr, int**mptr, double ***aptr, int *nmodesptr)
{
	FILE	*distparfile;
	int	*l, *m, nmodes, mode;
	char	line[1024];
	double	xorigin[2];
	double	**a;
	
	distparfile = fopen(filename, "r");
	if (!distparfile) 
		error_exit("getmodeamplitudes: failed to open distparfile\n");

	/* read the distortion parameters */
	while (fgets(line, 1024, distparfile)) {
		if (line[0] != '#')
				break;
	}
	if (1 != sscanf(line, "%d", &nmodes))
		error_exit("getmodeamplitudes: malformed distortion parameter file\n");
	if (!fgets(line, 1024, distparfile)) {
		error_exit("getmodeamplitudes: malformed distortion parameter file\n");
	}
	if (2 != sscanf(line, "%lf %lf", &(xorigin[0]), &(xorigin[1])))
		error_exit("getmodeamplitudes: malformed distortion parameter file\n");
	setorigin(xorigin);
	l = (int *) calloc(nmodes, sizeof(int));
	m = (int *) calloc(nmodes, sizeof(int));
	*aptr = a = (double **) calloc(2, sizeof(double *));
	a[0] = (double *) calloc(nmodes, sizeof(double));
	a[1] = (double *) calloc(nmodes, sizeof(double));
	mode = 0;
	while (mode < nmodes) {
		if (!fgets(line, 1024, distparfile))
			error_exit("getmodeamplitudes: malformed distortion parameter file\n");
		if (line[0] == '#')
			continue;
		if (4 != sscanf(line, "%d %d %lf %lf", l + mode, m + mode, a[0] + mode, a[1] + mode))		
			error_exit("getmodeamplitudes: malformed distortion parameter file\n");
		mode++;
	}
	*lptr = l;
	*mptr = m;
	*nmodesptr = nmodes;
}



void	setorigin(double *x)
{
	x0[0] = x[0];
	x0[1] = x[1];
}

void	modefunc_addargcomment(int argc, char *argv[])
{
	int	arg;

	argstring = (char *) calloc(1024, sizeof(char));
	for (arg = 0; arg < argc; arg++) {
		strcat(argstring, argv[arg]);
		strcat(argstring, " ");
	}	
}


int	getvars(FILE *lcpipe, int *nvar, char *vardef[], int *thesize)
{
	char	line[1024], word[10];
	int	contents = 0, gotvars = 0, pos, dpos, ndim, dim, i, itemsize, size = 0;

	*nvar = 0;
	while (!contents) {
		if (!fgets(line, 1024, lcpipe)) {
			return(0);
		}
		if (!strncmp(line, "# contents: ", 12)) {
			contents = 1;
		}
	}
	while (!gotvars) {
		if (!fgets(line, 1024, lcpipe)) {
			return(0);
		}
		if (!strncmp(line, "# text", 6)) {
			error_exit("getvars: catalogue must not include text object items!\n");
		}
		if (strncmp(line, "# number", 8)) {
			gotvars = 1;
		} else {
			vardef[*nvar] = (char *) calloc(128, sizeof(char));
			pos = 8;
			sscanf(line + pos, "%d%n", &ndim, &dpos);
			if (ndim >= 10) {
				return(0);
			}
			pos += dpos;
			sprintf(vardef[*nvar], "%d", ndim);
			itemsize = 1;
			for (i = 0; i < ndim; i++) {
				sscanf(line + pos, "%s%n", word, &dpos);
				pos += dpos;
				strcat(vardef[*nvar], " ");
				strcat(vardef[*nvar], word);
				sscanf(word, "%d", &dim);
				itemsize *= dim;
			}
			size += itemsize;
			sscanf(line + pos, "%s%n", word, &dpos);
			strcat(vardef[*nvar], " ");
			strcat(vardef[*nvar], word);

			(*nvar)++;
			if (*nvar == MODEFUNC_MAX_VARS) {
				return(0);
			}
		}
	}
	*thesize = size;
	return(1);
}