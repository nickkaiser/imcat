/*
 * powell_dbl.c
 *
 * original code from NR - converted to double and zero-base indexing
 * by Nick Kaiser
 *
 */


#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 200

void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret,
	double (*func)(double []))
{
	void linmin(double p[], double xi[], int n, double *fret,
		double (*func)(double []));
	int i,ibig,j;
	double del,fp,fptt,t,*pt,*ptt,*xit;

	pt=(double *) calloc(n, sizeof(double));
	ptt=(double *) calloc(n, sizeof(double));
	xit=(double *) calloc(n, sizeof(double));
	*fret=(*func)(p);
	for (j=0;j<n;j++) pt[j]=p[j];
	for (*iter=1;;++(*iter)) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			linmin(p,xit,n,fret,func);
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			free(xit);
			free(ptt);
			free(pt);
			return;
		}
		if (*iter == ITMAX) nrerror("powell exceeding maximum iterations.");
		for (j=0;j<n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,n,fret,func);
				for (j=0;j<n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
}
#undef ITMAX
#undef NRANSI
