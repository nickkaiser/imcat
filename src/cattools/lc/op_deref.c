/*
 * op_deref.c
 */

#include <stdio.h>
#include "../../catlib/cat.h"
#include "stack.h"
#include "error.h"

void	derefinit(item *operator)
{
	int	idim, index, dim0, ndim;
	item	*result, *op1, *op2;

	op2 = pop();
	op1 = pop();
	/* check the argument type */
	if ((op2->itype != NUM_TYPE) || (op2->ndim > 1)) {
		error_exit("derefinit: index format error\n");
	}
	index = (int) ((double *) (op2->addr))[0];
	dim0 = (op1->dim)[0];
	if (index < 0 || index >= dim0) {
		error_exit("deferinit: index out of range\n");
	}
	result = operator->next = (item *) calloc(1, sizeof(item));
	result->name = op1->name;
	result->itype = op1->itype;
	ndim = op1->ndim;
	if (ndim == 1) {
		result->ndim = 1;
		(result->dim)[0] = 1;
		switch (op1->itype) {
			case NUM_TYPE:
				result->addr = (void *) ((double *) op1->addr + index);
				break;
			case TEXT_TYPE:
				result->addr = (void *) ((char **) op1->addr + index);
				break;
			default:
				error_exit("derefinit: bad type\n");
				break;
		}
	} else {
		result->ndim = ndim - 1;
		for (idim = 0; idim < result->ndim; idim++) {
			(result->dim)[idim] = (op1->dim)[idim + 1];
		}
		result->addr = ((void **) (op1->addr))[index];
	}	
}



void	derefdoit(item *operator)
{
	item	*op1, *op2;

	op2 = pop();
	op1 = pop();
}


