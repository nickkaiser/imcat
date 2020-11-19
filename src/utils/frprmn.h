/*
 * frprmn.h
 *
 * copied from num rec + ansified by Nick Kaiser
 */

void frprmn(float p[], int n,float ftol,int *iter,float *fret,
float (*func)(float []), void (*dfunc)(float [], float []));
