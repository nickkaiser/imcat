/*
 * op_misc.c - if (a.k.a. ?) and enter
 */

#include <stdio.h>
#include "../../catlib/cat.h"
#include "stack.h"
#include "error.h"

void	ifinit(item *operator)
{
	item	*result, *op1, *op2, *op3;
	int	i;
	double	*cond;

	op3 = pop();
	if ((op3->itype != NUM_TYPE) || ((op3->dim)[0] != 1)) {
		error_exit("if: condition must be a single number\n");
	}
	cond = (double *) (op3->addr);
	op2 = pop();
	op1 = pop();
	if ((op1->itype != op2->itype) || (op1->ndim != 1) || (op2->ndim != 1)) {
		error_exit("if: 1st two arguments must be of same size, type\n");
	}
	for (i = 0; i < op1->ndim; i++) {
		if ((op1->dim)[i] != (op2->dim)[i]) {
			error_exit("if: 1st two arguments must be of same size, type\n");
		}
	}
	result = operator->next = copyitem(op1);
	allocitemcontents(result, &(result->addr), 0);
}

void	ifdoit(item *operator)
{
	item	*op1, *op2, *op3;
	double	*cond;

	op3 = pop();
	op2 = pop();
	op1 = pop();
	cond = (double *) (op3->addr);
	if (*cond != 0.0) {
		copyitemcontents(operator->next, (operator->next)->addr, op1, op1->addr, 0);
	} else {
		copyitemcontents(operator->next, (operator->next)->addr, op2, op2->addr, 0);
	}
}


void	enterinit(item *operator)
{
	item	*theop, *result;

	theop = pop();
	push(theop);
	result = operator->next = copyitem(theop);
	allocitemcontents(result, &(result->addr), 0);
}


void	enterdoit(item *operator)
{
	item	*theop;

	theop = pop();
	push(theop);
	copyitemcontents(operator->next, (operator->next)->addr, theop, theop->addr, 0);
}
