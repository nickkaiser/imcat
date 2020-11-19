/*
 * stack.c - rpn stack functions for item arithmetic
 *
 */

#include <stdio.h>
#include "../../catlib/cat.h"
#include "stack.h"
#include "error.h"

#define STACK_DEPTH 100

item	*itemstack[STACK_DEPTH];
int	stackpos = 0;

void	push(item *theitem)
{
/*
	printf("pushed: ");
	writename(theitem);
	writeitem(theitem, theitem->addr, 0);
	printf("\n");
*/
	if (stackpos < STACK_DEPTH) {
		itemstack[stackpos++] = theitem;
	} else {
		error_exit("push: stack full\n");
	}
}


item	*pop(void)
{
/*
	printf("popped: ");
	writename(itemstack[stackpos - 1]);
	writeitem(itemstack[stackpos - 1], itemstack[stackpos - 1]->addr, 0);
	printf("\n");
*/
	if (stackpos > 0) {
		return(itemstack[--stackpos]);
	} else {
		error_exit("pop: stack emtpy\n");
	}
}
