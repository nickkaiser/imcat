/*
 * getop.c
 *
 */

#include <stdio.h>
#include <limits.h>
#include "getop.h"
#include "../../utils/error.h"
#include "../../imlib/fits.h"

#define MAGIC	FLOAT_MAGIC
extern int 	nim;
static int	N1;


static char     *num0funcname[N_NUM0_OPS] = {
"rand", "x", "y", "xp", "yp", "if", "?", "enter", "grand"};

static char     *num1funcname[N_NUM1_OPS] = {
"acos", "asin", "atan", "ceil", "cos", "cosh", "exp", "fabs", "floor",
"log", "log10", "sin", "sinh", "sqrt", "tan", "tanh", "!", "j0", "j1", "y0", "y1", "swap"
};

static char     *num2funcname[N_NUM2_OPS] = {
"mult", "*",    "+",  "/",   "-", "atan2", "pow", "fmod",
">", ">=", "<", "<=", "!=", "==", "jn", "yn", "max", "min"
};


void	set_getop_N1(int n1)
{
	N1 = n1;
}


op	*getop(char **sptr)
{
	char	*theword, *wordbase;
	int	nchars, fitsnum;
	double	number;
	op	*theop;
	int	ix;
	
	wordbase = theword = getword(sptr);

	if (!theword) {
		return (NULL);
	}
	theop = (op *) calloc(1, sizeof(op));
	theop->data = (double *) calloc(N1, sizeof(double));

	/* check for MAGIC */
	if (!strcmp(theword, "MAGIC")) {
		theop->type = CONSTANT_TYPE;
		theop->opno = 0;
		for (ix = 0; ix < N1; ix ++) {
			theop->data[ix] = MAGIC;
		}
		return (theop);
	}
	/* check for a constant number */
	/* used to "%lg %n" here */
	if (1 == sscanf(theword, "%lg%n", &number, &nchars)) {
		if (nchars != strlen(theword)) {
			fprintf(stderr, "getop: can't understand %s\n", theword);
			exit(-1);
		} else {
			theop->type = CONSTANT_TYPE;
			theop->opno = 0;
			for (ix = 0; ix < N1; ix ++) {
				theop->data[ix] = number;
			}
			return (theop);
		}
	}
	/* now check for variables, operators */
	switch (theword[0]) {
		case '%':			/* fits file */
			theword++;
			if (!strlen(theword)) {
				error_exit("getop: expression format error\n");
			}
			/* used to have "%d %n" here */
			if (1 == sscanf(theword, "%d%n", &fitsnum, &nchars)) {
				if (nchars != strlen(theword)) {
					fprintf(stderr, "getop: can't understand %s\n", wordbase);
					exit(-1);
				}
				if ((fitsnum < 1) || (fitsnum > nim)) {
					error_exit("getop: fits image out of bounds\n");
				}
				theop->type = IM_VALUE_TYPE;	
				theop->opno = fitsnum - 1;
				return (theop);
			}
			break;
		default:
			return(newop(theword));
			break;		
	}
	fprintf(stderr, "getop: can't understand %s\n", wordbase);
	exit(-1);
}

/* basically want to return tolens separated by white space but slightly */
/* tricky as we want to convert "[y]" to "y deref" */
char	*getword(char **sptr)
{
	char	c, *string, *theword;
	int	start, end, started, ended, len;

	string = *sptr;
	start = started = 0;
	/* skip over any leading space */
	while(!started) {
		c = string[start];
		switch (c) {
			case '\0':
				return (NULL);
				break;
			case ' ':
			case '\t':
			case '[':
				start++;
				break;
			default:
				started = 1;
				break;
		}
	}
	/* find the char after the end */
	end = start;
	ended = 0;
	while (!ended) {
		c = string[end];
		switch (c) {
			case '\0':
			case  ' ':
			case '\t':
			case  '[':
				ended = 1;
				break;
			default:
				end++;
				break;
		}
	}
	if ((end != (start+1)) && (string[end-1] == ']')) {
		end--;
	}
	*sptr += end;
	if (end == start)
		return (NULL);
	len = end - start;
	if ((len == 1) && (string[start] == ']')) {
		theword = (char *) calloc(6, sizeof(char));
		strcpy(theword,  "deref");
	} else {
		theword = (char *) calloc(len + 1, sizeof(char));
		strncpy(theword,  string + start, len);
	}
	return (theword);	
}



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




