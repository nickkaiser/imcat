/*
 * operators.c
 *
 */

#include <stdio.h>
#include <math.h>
#include "stack.h"
#include "getop.h"
#include "operators.h"
#include "op_math.h"
#include "op_logic.h"
#include "ic.h"
/*
#include "op_rand.h"
#include "op_deref.h"
#include "op_vector.h"
#include "op_dot.h"
#include "error.h"
*/

/* nonadic numeric functions */
#define N_NUM0_OPS	8
static char 	*num0funcname[N_NUM0_OPS] = {
"rand", "x", "y", "xp", "yp", "if", "?", "enter"};
double	(*num0funcarray[N_NUM0_OPS])(void) = {
drand48, x, y, xp, yp, If, If, enter};

/* monadic numeric functions */
#define N_NUM1_OPS	17
static char 	*num1funcname[N_NUM1_OPS] = {
"acos", "asin", "atan", "ceil", "cos", "cosh", "exp", "fabs", "floor",
"log", "log10", "sin", "sinh", "sqrt", "tan", "tanh", "!"
};
double	(*num1funcarray[N_NUM1_OPS])(double) = {
 acos,   asin,   atan,   ceil,   cos,   cosh,   exp,   fabs,   floor,
  log,  log10,   sin,  sinh,  sqrt,  tan,  tanh, not
};

/* dyadic numeric functions */
#define N_NUM2_OPS	13
static char 	*num2funcname[N_NUM2_OPS] = {
"mult", "*",    "+",  "/",   "-", "atan2", "pow", "fmod",
">", ">=", "<", "<=", "!="
};
double	(*num2funcarray[N_NUM2_OPS])(double, double) = {
times, times, plus, divide, minus, atan2,   pow,   fmod,
gt, ge, lt, le, ne
};




op		*newop(char *thename)
{
	int	iop;
	op	*theop;

	/* first handle the numeric operators */
	for (iop = 0; iop < N_NUM2_OPS; iop++) {
		if (!strcmp(num2funcname[iop], thename)) {
			break;
		}
	}
	if (iop < N_NUM2_OPS) {
		theop = (op *) calloc(1, sizeof(op));
		theop->type = NUM2_FUNC_TYPE;
		theop->opno = iop;
		return (theop);
	}

	for (iop = 0; iop < N_NUM1_OPS; iop++) {
		if (!strcmp(num1funcname[iop], thename)) {
			break;
		}
	}
	if (iop < N_NUM1_OPS) {
		theop = (op *) calloc(1, sizeof(op));
		theop->type = NUM1_FUNC_TYPE;
		theop->opno = iop;
		return (theop);
	}

	for (iop = 0; iop < N_NUM0_OPS; iop++) {
		if (!strcmp(num0funcname[iop], thename)) {
			break;
		}
	}
	if (iop < N_NUM0_OPS) {
		theop = (op *) calloc(1, sizeof(op));
		theop->type = NUM0_FUNC_TYPE;
		theop->opno = iop;
		return (theop);
	}

	fprintf(stderr, "newopitem: unknown function: %s\n", thename);
	exit(0);
}







