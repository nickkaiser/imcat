/*
 * ipbuff.h
 */

double	**readdoublebuff(int width, FILE *stream, int *np);
/* returns address of array of *np pointers to input records of size width * sizeof(double) */
