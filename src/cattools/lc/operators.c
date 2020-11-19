/*
 * operators.c
 *
 */

#include <stdio.h>
#include <math.h>
#include "../../catlib/cat.h"
#include "stack.h"
#include "operators.h"
#include "op_math.h"
#include "op_logic.h"
#include "op_deref.h"
#include "op_vector.h"
#include "op_dot.h"
#include "op_rand.h"
#include "op_misc.h"
#include "op_lintrans.h"
#include "error.h"
#include "op_strings.h"

/* monadic numeric functions */
#define N_NUM1_OPS	21
static char 	*num1funcname[N_NUM1_OPS] = {
"acos", "asin", "atan", "ceil", "cos", "cosh", "exp", "fabs", "floor",
"log", "log10", "sin", "sinh", "sqrt", "tan", "tanh", "!", "j0", "j1", "y0", "y1"
};
double	(*num1funcarray[N_NUM1_OPS])(double) = {
 acos,   asin,   atan,   ceil,   cos,   cosh,   exp,   fabs,   floor,
  log,  log10,   sin,  sinh,  sqrt,  tan,  tanh, not, j0, j1, y0, y1
};

/* dyadic numeric functions */
#define N_NUM2_OPS	18
static char 	*num2funcname[N_NUM2_OPS] = {
"mult",  "*",    "+",  "/",   "-", "atan2", "pow", "fmod",
">", ">=", "<", "<=", "==", "!=", "and", "or", "max", "min"
};
double	(*num2funcarray[N_NUM2_OPS])(double, double) = {
times, times, plus, divide, minus, atan2,   pow,   fmod,
gt, ge, lt, le, eq, ne, and, or, max, min
};


/* generic functions */
#define N_GEN_OPS	18
static char *genfuncname[N_GEN_OPS] = {
"deref", "vector", "dot", "rand", "grand", "if", "?", "enter", "lintrans", "vsub", "vadd", "vscale", "vshift", 
	"inverse", "msub", "madd", "mscale", "eq"};
void	(*genfuncinit[N_GEN_OPS])(item *) = {
derefinit, vectorinit, dotinit, randinit, grandinit, ifinit, ifinit, enterinit, lintransinit, vsubinit, vaddinit, vscaleinit, vshiftinit, 
	inverseinit, msubinit, maddinit, mscaleinit, streqinit};
void	(*genfuncdoit[N_GEN_OPS])(item *) = {
derefdoit, vectordoit, dotdoit, randdoit, granddoit, ifdoit, ifdoit, enterdoit, lintransdoit, vsubdoit, vadddoit, vscaledoit, vshiftdoit, 
	inversedoit, msubdoit, madddoit, mscaledoit, streqdoit};

item		*newopitem(char *thename)
{
	int	iop;
	item	*theitem;

	/* first handle the numeric operators */
	for (iop = 0; iop < N_NUM2_OPS; iop++) {
		if (!strcmp(num2funcname[iop], thename)) {
			break;
		}
	}
	if (iop < N_NUM2_OPS) {
		theitem = newitem(thename, NUM2_OP_TYPE, 1, 1);
		(theitem->dim)[0] = iop;
		theitem->addr = (void *) num2funcarray[iop];
		return (theitem);
	}
	for (iop = 0; iop < N_NUM1_OPS; iop++) {
		if (!strcmp(num1funcname[iop], thename)) {
			break;
		}
	}
	if (iop < N_NUM1_OPS) {
		theitem = newitem(thename, NUM1_OP_TYPE, 1, 1);
		(theitem->dim)[0] = iop;
		theitem->addr = (void *) num1funcarray[iop];
		return (theitem);
	}
	for (iop = 0; iop < N_GEN_OPS; iop++) {
		if (!strcmp(genfuncname[iop], thename)) {
			break;
		}
	}
	if (iop < N_GEN_OPS) {
		theitem = newitem(thename, GEN_OP_TYPE, 1, 1);
		(theitem->dim)[0] = iop;
		return (theitem);
	}
	fprintf(stderr, "newopitem: unknown function: %s\n", thename);
	exit(-1);
}


void	num1func(item *operand, item *operator)
{
	double	*input, *output;
	
	/* check that we have the right kind of operand */
	if ((operand->itype != NUM_TYPE) || (operand->ndim != 1) || ((operand->dim)[0] != 1)) {
		error_exit("num1func: argument type error\n");
	}

	/* allocate space for result if it doesn't already exist */
	if (!(operator->next)) {
		operator->next = newitem("TEMP", NUM_TYPE, 1, 1);
		allocitemcontents(operator->next, &((operator->next)->addr), 0);
	}

	input = (double *) (operand->addr);
	output = (double *) ((operator->next)->addr);
	output[0] = num1funcarray[(operator->dim)[0]](input[0]);
}

void	num2func(item *operand1, item *operand2, item *operator)
{
	double	*input1, *input2, *output;

	/* check that we have the right kind of operands */
	if ((operand1->itype != NUM_TYPE) || (operand1->ndim != 1) || ((operand1->dim)[0] != 1) ||
		(operand2->itype != NUM_TYPE) || (operand2->ndim != 1) || ((operand2->dim)[0] != 1)) {
		error_exit("num2func: argument type error\n");
	}

	/* allocate space for result if it doesn't already exist */
	if (!(operator->next)) {
		operator->next = newitem("TEMP", NUM_TYPE, 1, 1);
		allocitemcontents(operator->next, &((operator->next)->addr), 0);
	}

	input1 = (double *) (operand1->addr);
	input2 = (double *) (operand2->addr);
	output = (double *) ((operator->next)->addr);
	output[0] = num2funcarray[(operator->dim)[0]](input1[0], input2[0]);
}




/* generic function definition */
void	genfunc(item *operator)
{
	int	idim, index;

	/* first time around we need to */
	/* check that we have the right kind of operand, figure type of result and */
	/* allocate space for result */
	if (!(operator->next)) {
		genfuncinit[(operator->dim)[0]](operator);		
	} else {
		genfuncdoit[(operator->dim)[0]](operator);		
	}
}







