/*
 * amoeba_dbl.c
 *
 * original code from NR converted to double and zero-based indexing by Nick Kaiser
 */


#include <math.h>
#define NRANSI
#include "nrutil.h"
#define NMAX 500000
#define GET_PSUM \
					for (j=0;j<ndim;j++) {\
					for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j];\
					psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void amoeba(double **p, double y[], int ndim, double ftol,
	double (*funk)(double []), int *nfunk)
{
	double amotry(double **p, double y[], double psum[], int ndim,
		double (*funk)(double []), int ihi, double fac);
	int i,ihi,ilo,inhi,j,mpts=ndim+1;
	double rtol,sum,swap,ysave,ytry,*psum;

	psum=(double *) calloc(ndim, sizeof(double));
	*nfunk=0;
	GET_PSUM
	for (;;) {
		ilo=1;
		ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
		for (i=0;i<mpts;i++) {
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi) inhi=i;
		}
		rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
		if (rtol < ftol) {
			SWAP(y[1],y[ilo])
			for (i=0;i<ndim;i++) SWAP(p[1][i],p[ilo][i])
			break;
		}
		if (*nfunk >= NMAX) nrerror("NMAX exceeded");
		*nfunk += 2;
		ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
		if (ytry <= y[ilo])
			ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);
		else if (ytry >= y[inhi]) {
			ysave=y[ihi];
			ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
			if (ytry >= ysave) {
				for (i=0;i<mpts;i++) {
					if (i != ilo) {
						for (j=0;j<ndim;j++)
							p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
						y[i]=(*funk)(psum);
					}
				}
				*nfunk += ndim;
				GET_PSUM
			}
		} else --(*nfunk);
	}
	free(psum);
}
#undef SWAP
#undef GET_PSUM
#undef NMAX
#undef NRANSI
