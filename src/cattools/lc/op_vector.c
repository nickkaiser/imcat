/*
 * op_vector.c
 */

#include <stdio.h>
#include "../../catlib/cat.h"
#include "stack.h"
#include "error.h"



void	vectorinit(item *operator)
{
	item	*result, *op;
	int	itype, ndim, *dim, i, ii, idim;
	void	**addr;
	int	ncomp, dim0;

	/* get the number of components */
	op = pop();
	if ((op->itype != NUM_TYPE) || (op->ndim != 1) || ((op->dim)[0] != 1)) {
		error_exit("vectorinit: non-numeric number of components\n");
	}
	ncomp = (int) ((double *) (op->addr))[0];

	/* get the first argument */
	op = pop();
	itype = op->itype;
	ndim = op->ndim;
	dim = op->dim;
	dim0 = (op->dim)[0];
	if (dim0 > 1) {
		addr = (void **) calloc(ncomp, sizeof(void *));
	} else {
		switch (itype) {
			case NUM_TYPE:
				addr = (void **) calloc(ncomp, sizeof(double));
				break;
			case TEXT_TYPE:
				addr = (void **) calloc(ncomp, sizeof(char *));
				break;
			default:
				error_exit("vectorinit: bad type\n");
				break;
		}
	}
	for (i = 0; i < ncomp; i++) {
		ii = ncomp - i - 1;
		if (i)
			op = pop();
		if (dim0 > 1) {
			addr[ii] = op->addr;
		} else {
			switch (itype) {
				case NUM_TYPE:
					((double *) addr)[ii] = *((double *) (op->addr));
					break;
				case TEXT_TYPE:
					((char **) addr)[ii] = *((char **) (op->addr));
					break;
				default:
					error_exit("vectorinit: bad type\n");
					break;
			}
		}
		if ((op->itype != itype) || (op->ndim != ndim)) {
			error_exit("vectorinit: non identical component types\n");
		}
		for (idim = 0; idim < ndim; idim++) {
			if (dim[idim] != (op->dim)[idim]) {
				error_exit("vectorinit: non identical component types\n");
			}
		}
	}	
	result = operator->next = (item *) calloc(1, sizeof(item));
	result->name = op->name;
	result->itype = op->itype;
	if (dim0 > 1) {
		result->ndim = 1 + op->ndim;
		if (result->ndim >= MAX_DIMS) {
			error_exit("vectorinit: too many dimensions\n");
		}
	} else {
		result->ndim = op->ndim;
	}
	(result->dim)[0] = ncomp;
	if (dim0 > 1) {
		for (idim = 0; idim < ndim; idim++) {
			(result->dim)[idim + 1] = (op->dim)[idim];
		}
	}
	result->addr = (void *) addr;
}


void	vectordoit(item *operator)
{
	int	i, ii, ncomp, ndim;
	item	*op;
	void	**addr;

	ndim = (operator->next)->ndim;
	ncomp = ((operator->next)->dim)[0];
	addr = (void **) ((operator->next)->addr);

	pop();
	for (i = 0; i < ncomp; i++) {
		ii = ncomp - i - 1;
		op = pop();
		if (ndim == 1) {
			switch (op->itype) {
				case NUM_TYPE:
					((double *) addr)[ii] = *((double *) (op->addr));
					break;
				case TEXT_TYPE:
					((char **) addr)[ii] = *((char **) (op->addr));
					break;
				default:
					error_exit("vectorinit: bad type\n");
					break;
			}
		}
	}
}
