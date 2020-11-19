/*
 * lmodelcalculus.c - functions for differentiating, integrating lmodels
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "catlib/cat.h"
#include "utils/error.h"
#include "utils/args.h"
#include "lmodel.h"


/* create derivative of an lmodel */
lmodel	*difflmodel(lmodel *ipmodel, int integratemode)
{
	lmodel	*opmodel;
	item	*aitem, the_aitem;
	int	i, j, mi, mo, equal, kzero;
	int	lmin, lmax, xdim;
	double	pp, kk;

	/* local copy of xdim */
	xdim = ipmodel->xdim;

	/* copy the input a-item */
	the_aitem = *(ipmodel->aitem);
	aitem = &the_aitem;
	/* and vectorise or devectorize it */
	if (integratemode) {
		(aitem->ndim)--;
	} else {
		aitem->dim[aitem->ndim] = xdim; 
		(aitem->ndim)++;
	}

	/* generate the output model opmodel */
	switch (ipmodel->modeltype) {
		case POLYNOMIAL_LMODEL:
			if (integratemode) {
				lmin = ipmodel->lmin + 1;
				lmax = ipmodel->lmax + 1;
			} else {
				lmax = ipmodel->lmax - 1;
				lmin = (ipmodel->lmin ? ipmodel->lmin - 1 : 0);
			}
			opmodel = newpolylmodel(ipmodel->xitem, aitem, lmin, lmax);
			break;
		case ZERNIKE_LMODEL:
			error_exit("difflmodel: can't do Zernike models\n");
			break;
		case FOURIER_LMODEL:
			opmodel = newfourierlmodel(ipmodel->xitem, aitem, ipmodel->kmin, ipmodel->kmax, ipmodel->boxlimits, ipmodel->lbox);
			break;
		default:
			error_exit("difflmodel: unrecognized model type\n");
	}

	flatten_a(ipmodel);
	flatten_a(opmodel);

	for (mi = 0; mi < ipmodel->nmodes; mi++) {
	    for (j = 0; j < xdim; j++) {
		switch (ipmodel->modeltype) {
			case POLYNOMIAL_LMODEL:
				for (i = 0; i < xdim; i++) {
					opmodel->pp[i] = (int) ipmodel->p[mi][i];
				}
				if (integratemode) {
					opmodel->pp[j]++;
					pp = 0.0;
					for (i = 0; i < xdim; i++) {
						pp += opmodel->pp[i] * opmodel->pp[i]; 
					}
				} else {
					opmodel->pp[j]--;
				}
				for (mo = 0; mo < opmodel->nmodes; mo++) {
					equal = 1;
					for (i = 0; i < xdim; i++) {
						if (opmodel->p[mo][i] != (double) (opmodel->pp[i])) {
							equal = 0;
							break;
						}
					}
					if (equal) {
						for (i = 0; i < (integratemode ? opmodel->asize : ipmodel->asize); i++) {
							if (integratemode) {
								if (pp > 0) {
									opmodel->aflat[i][mo] += ipmodel->aflat[i * xdim + j][mi] / pp;
								}
							} else {
								opmodel->aflat[i * xdim + j][mo] = ipmodel->p[mi][j] * ipmodel->aflat[i][mi];
							}
						}
					}
				}
				break;
			case ZERNIKE_LMODEL:
				error_exit("difflmodel: can't do Zernike models\n");
				break;
			case FOURIER_LMODEL:
				kzero = 1;
				kk = 0.0;
				for (i = 0; i < xdim; i++) {
					if (ipmodel->k[mi][i] != 0.0) {
						kzero = 0;
					}
					kk += ipmodel->k[mi][i] * ipmodel->k[mi][i];
				}
				for (i = 0; i < (integratemode ? opmodel->asize : ipmodel->asize); i++) {
					if (kzero) {
						opmodel->aflat[i][mi] = 0.0;
					} else {
						mo = (ipmodel->i[mi] != 0.0 ? mi - 1 : mi + 1);
						if (integratemode) {
							opmodel->aflat[i][mo] += ipmodel->lbox * ipmodel->k[mi][j] * ipmodel->aflat[i * xdim + j][mi] / (2 * M_PI * kk);
						} else {
							opmodel->aflat[i * xdim + j][mo] = 2 * M_PI * ipmodel->k[mi][j] * ipmodel->aflat[i][mi] / ipmodel->lbox;
						}
					}
				}
				break;
			default:
				error_exit("difflmodel: unrecognized model type\n");
		}
	    }
	}

	unflatten_a(opmodel);
	return(opmodel);
}