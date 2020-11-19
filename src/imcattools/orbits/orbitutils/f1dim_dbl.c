/*
 * f1dim_dbl.c
 *
 * original code from NR - converted to double and zero-base indexing
 * by Nick Kaiser
 *
 */


#define NRANSI
#include "nrutil.h"

extern int ncom;
extern double *pcom,*xicom,(*nrfunc)(double []);

double f1dim(double x)
{
	int j;
	double f,*xt;

	xt=(double *) calloc(ncom, sizeof(double));
	for (j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free(xt);
	return f;
}
#undef NRANSI
