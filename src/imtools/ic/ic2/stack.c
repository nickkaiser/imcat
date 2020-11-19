/*
 * stack.c - rpn stack functions for float arithmetic
 *
 */

#include <stdio.h>
#include "stack.h"
#include "error.h"

#define STACK_DEPTH 100

double	stack[STACK_DEPTH];
int	stackpos = 0;

void	push(double number)
{
/*
	printf("debug: pushed %lf\n", number);
*/
	if (stackpos < STACK_DEPTH) {
		stack[stackpos++] = number;
	} else {
		error_exit("push: stack full\n");
	}
}


double	pop(void)
{
/*
	printf("debug: popping %lf\n", stack[stackpos - 1]);
*/
	if (stackpos > 0) {
		return(stack[--stackpos]);
	} else {
		error_exit("pop: stack emtpy\n");
	}
}


void	resetstack(void)
{
	stackpos = 0;
}
