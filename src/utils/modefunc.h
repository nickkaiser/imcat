/*
 * modefunc.h
 */

#define MODEFUNC_MAX_VARS 100

double	f(int l, int m, double *x);
void	get2Dpolymodel(char *filename, int **lptr, int**mptr, int *asize, double ***aptr, int *nmodesptr, int *nvar, char *vardef[], char **xname);
void	getmodeamplitudes_txt(char *filename, int **lptr, int**mptr, double ***aptr, int *nmodesptr);
void	getmodeamplitudes_lc(char *filename, int **lptr, int**mptr, int *asize, double ***aptr, int *nmodesptr, int *nvar, char *vardef[], char **xname);
void	setorigin(double *x);
int	write2Dpolymodel(char *parfile, int nmodes, int *l, int *m, int asize, double **a, int nvar, char *vardef[], char *xvar);
void	modefunc_addargcomment(int argc, char *argv[]);
int	getvars(FILE *lcpipe, int *nvar, char *vardef[], int *size);

