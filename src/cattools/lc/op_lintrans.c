/*
 * op_lintrans.c - linear transformation: as in '%x phi00 phi01 phi10 phi11 lintrans'
 */

#include <stdio.h>
#include "../../catlib/cat.h"
#include "stack.h"
#include "error.h"

void	lintransinit(item *operator)
{
	item	*result, *op[4], *vec;
	int	i, m, n;
	double	*xin, *xout, *phi[2][2];

	for (i = 3; i >= 0; i--) {
		op[i] = pop();
		if ((op[i]->itype != NUM_TYPE) || (op[i]->ndim != 1) || ((op[i]->dim)[0] != 1)) {
			error_exit("lintrans: bad argument\n");
		}
	}
	vec = pop();
	if ((vec->itype != NUM_TYPE) || (vec->ndim != 1) || ((vec->dim)[0] != 2)) {
		error_exit("lintrans: bad argument\n");
	}
	phi[0][0] = (double *) op[0]->addr;
	phi[0][1] = (double *) op[1]->addr;
	phi[1][0] = (double *) op[2]->addr;
	phi[1][1] = (double *) op[3]->addr;
	xin = (double *) vec->addr;
	result = operator->next = newitem("LINTRANS_TEMP", NUM_TYPE, 1, 2);
	xout = (double *) calloc(2, sizeof(double));
	result->addr = (void *) xout;
	for (m = 0; m < 2; m++) {
		xout[m] = 0;
		for (n = 0; n < 2; n++) {
			xout[m] += *(phi[m][n]) * xin[n];
		}
	}
}

void	lintransdoit(item *operator)
{
	item	*result, *op[4], *vec;
	int	i, m, n;
	double	*xin, *xout, *phi[2][2];

	for (i = 3; i >= 0; i--) {
		op[i] = pop();
	}
	vec = pop();
	phi[0][0] = (double *) op[0]->addr;
	phi[0][1] = (double *) op[1]->addr;
	phi[1][0] = (double *) op[2]->addr;
	phi[1][1] = (double *) op[3]->addr;
	xin = (double *) vec->addr;
	result = operator->next;
	xout = (double *) result->addr;
	for (m = 0; m < 2; m++) {
		xout[m] = 0;
		for (n = 0; n < 2; n++) {
			xout[m] += *(phi[m][n]) * xin[n];
		}
	}
}
