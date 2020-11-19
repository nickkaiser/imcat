/*
 * rpn.c
 */


#include <stdio.h>
#include "../../catlib/cat.h"
#include "getop.h"
#include "stack.h"
#include "operators.h"
#include "rpn.h"
#include "error.h"

#define MAX_OPS	1024

rpnfunc		*newrpnfunction(char *name, char *string)
{
	item	*theop, *op2, **op, *result;
	int	nop = 0;
	char	*temp;
	rpnfunc *therpnfunc;

	therpnfunc = (rpnfunc *) calloc(1, sizeof(rpnfunc));
	if (name) {
		therpnfunc->name = (char *) calloc(1 + strlen(name), sizeof(char));
		strcpy(therpnfunc->name, name);
	}	

	op = (item **) calloc(MAX_OPS, sizeof(item));
	therpnfunc->oplist = op;

	/* get the operands and operators into an array of items */
	while (theop = getop(&string)) {
			op[nop++] = theop;
/*
			writename(theop);
			printf("address: %ld\n", (long) (theop->addr));
			writeitem(theop, theop->addr, 0);
			printf("\n");
*/
	}
	therpnfunc->nop = nop;
	
	return (therpnfunc);
}


void	evalrpnfunction(rpnfunc *therpnfunc)
{
	int	iop, nop;
	item	**op, *theop, *result, *op2;


	nop = therpnfunc->nop;
	op =  therpnfunc->oplist;

	/* feed the items through the stack machinery */
	for (iop = 0; iop < nop; iop++) {
		theop = op[iop];
		switch (theop->itype) {
			case NUM_TYPE:
			case TEXT_TYPE:
				push(theop);
				break;
			case NUM1_OP_TYPE:
				num1func(pop(), theop);
				result = theop->next;
				push(result);
				break;
			case NUM2_OP_TYPE:
				op2 = pop();
				num2func(pop(), op2, theop);
				result = theop->next;
				push(result);
				break;
			case GEN_OP_TYPE:
				genfunc(theop);
				result = theop->next;
				push(result);
				break;				
			default:
				break;
		}
		theop = theop->next;
	}
	result = pop();
	if (!(therpnfunc->result)) {
		therpnfunc->result = result;
		(therpnfunc->result)->name = therpnfunc->name;
	}
}
