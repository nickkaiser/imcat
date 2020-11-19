/*
 * vectors.c - definitions of some functions which operate on vectors
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"

void	assign(double *vec, double v0, double v1, double v2)
{
	vec[0] = v0;
	vec[1] = v1;
	vec[2] = v2;
}

double	length(double *vec)
{
	int	i;
	double	result;

	result = 0.0;
	for (i = 0; i < 3; i++) {
		result += vec[i] * vec[i];
	}
	return (sqrt(result));
}

void	diff(double *dst, double *src1, double *src2)
{
	int	i;

	for (i = 0; i < 3; i++) {
		dst[i] = src1[i] - src2[i];
	}
}

void	add(double *dst, double *src1, double *src2)
{
	int	i;

	for (i = 0; i < 3; i++) {
		dst[i] = src1[i] + src2[i];
	}
}

void	scale(double *vec, double scalefactor)
{
	int	i;

	for (i = 0; i < 3; i++) {
		vec[i] *= scalefactor;
	}
}

void	copy(double *dst, double *src)
{
	int	i;

	for (i = 0; i < 3; i++) {
		dst[i] = src[i];
	}
}

void	printvec(double *vec, char *name)
{
	int	i;

	fprintf(stdout, "%10s", name);
	for (i = 0; i < 3; i++) {
		fprintf(stdout, " %14.10lf", vec[i]);
	}
	fprintf(stdout, "\n");
}

void	fprintvec(double *vec, FILE *stream)
{
	int	i;

	for (i = 0; i < 3; i++) {
		fprintf(stream, " %14.10lf ", vec[i]);
	}
}

double	dot(double *vec1, double *vec2)
{
	int	i;
	double	result;

	result = 0.0;
	for (i = 0; i < 3; i++) {
		result += vec1[i] * vec2[i];
	}
	return (result);
}

void	getperp(double *dst, double *src, double *n)
{
	int	i;
	double	srcdotn;

	srcdotn = dot(src, n);
	for (i = 0; i < 3; i++) {
		dst[i] = src[i] - srcdotn * n[i];
	}
}

void	rotx(double *vec, double angle)
{
	/* rotate anticlockwise about x axis */
	double	x, y, z;

	x = vec[0];
	y = vec[1];
	z = vec[2];
	
	assign(vec, x, y * cos(angle) - z * sin(angle), y * sin(angle) + z * cos(angle));
}

void	roty(double *vec, double angle)
{
	/* rotate anticlockwise about y axis */
	double	x, y, z;

	x = vec[0];
	y = vec[1];
	z = vec[2];
	
	assign(vec, x * cos(angle) + z * sin(angle), y, z * cos(angle) - x * sin(angle));
}

void	rotz(double *vec, double angle)
{
	/* rotate anticlockwise about z axis */
	double	x, y, z;

	x = vec[0];
	y = vec[1];
	z = vec[2];
	
	assign(vec, x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle), z);
}

void	cross(double *dst, double *a, double *b)
{
	dst[0] = a[1] * b[2] - b[1] * a[2];
	dst[1] = a[2] * b[0] - b[2] * a[0];
	dst[2] = a[0] * b[1] - b[0] * a[1];
}
