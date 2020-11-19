/*
 * radecio.c -- functions to convert hms <=> decimal etc.
 */

#include <stdio.h>
#include <math.h>
#include "utils/error.h"
#include "radecio.h"

int	decimaltoxms(double angle, char *hmsstring)
{
	double 	h, m, s;
	int 	sign;

	sign = (angle < 0.0 ? -1 : 1);
	angle *= sign;
	angle = 60.0 * modf(angle, &h);
	s = 60.0 * modf(angle, &m);
	if (sign > 0) {
		sprintf(hmsstring, "%d:%d:%2.4lf", (int) h, (int) m, s);
	} else {
		sprintf(hmsstring, "-%d:%d:%2.4lf", (int) h, (int) m, s);
	}	
}


int	xmstodecimal(char *hmsstring, double *result)
{
	int	h, m, sign;
	double	s;

	if (hmsstring[0] == '-') {
		sign = -1;
		hmsstring = hmsstring + 1;
	} else {
		sign = 1;
	}
	if (3 != sscanf(hmsstring, "%d:%d:%lf", &h, &m, &s)) {
		return(0);
	} 
	if (m < 0 || m > 59 || s < 0.0 || s > 60.0) {
		error_exit("xmstodecimal: illegal minutes or seconds\n");
	}
	*result = sign * (h + m / 60.0 + s / 3600.0);
	return(1);
}


double	getangle(char *argstring, int *type)
{
	double	angle;
	int	pos;

	*type = DECIMAL_ANGLE_TYPE;
	if (1 != sscanf(argstring, "%lf%n", &angle, &pos)) {
		error_exit("getangle: can't decipher string\n");
	}

	if (pos != strlen(argstring)) {
		*type = HMS_ANGLE_TYPE;
		if (!xmstodecimal(argstring, &angle)) {
			error_exit("getangle: can't decipher string\n");
		}
	}

	return (angle);
}

double	getra(char *argstring, int *type)
{
	double	ra;

	ra = getangle(argstring, type);
	if (*type == HMS_ANGLE_TYPE) {
		ra *= 15.0;
	}
	if (ra > 360.0 || ra <= -360.0) {
		error_exit("getra: ra out of range -360.0 - 360.0\n");
	}

	return (ra);
}

double	getdec(char *argstring, int *type)
{
	double	dec;

	dec = getangle(argstring, type);
	if (dec > 90.0 || dec < -90.0) {
		error_exit("getdec: dec out of range -90.0 - 90.0\n");
	}

	return (dec);
}


