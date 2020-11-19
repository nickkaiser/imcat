/*
 * f1dim.c
 *
 * copied from num rec + ansified by Nick Kaiser
 */

#include <stdio.h>
#include "f1dim.h"
#include "nrutil.h"

extern int ncom;	/* defined in LINMIN */
extern float *pcom,*xicom,(*nrfunc)();

float f1dim(float x)
{
	int j;
	float f,*xt,*vector();

	xt=vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_vector(xt,1,ncom);
	return f;
}

