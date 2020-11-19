/*
 * radecio.h -- functions to convert hms <=> decimal
 */

#define	HMS_ANGLE_TYPE 0
#define DECIMAL_ANGLE_TYPE 1

int	decimaltoxms(double angle, char *hmsstring);
int	xmstodecimal(char *hmsstring, double *result);
double	getangle(char *argstring, int *type);
double	getra(char *argstring, int *type);
double	getdec(char *argstring, int *type);

