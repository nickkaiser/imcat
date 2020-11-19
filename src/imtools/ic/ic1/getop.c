/*
 * getop.c
 *
 */

#include <stdio.h>
#include <limits.h>
#include "getop.h"
#include "operators.h"
#include "../../utils/error.h"
#include "../magic.h"

extern int nim;


op	*getop(char **sptr)
{
	char	*theword, *wordbase;
	int	nchars, fitsnum;
	double	number;
	op	*theop;
	
	wordbase = theword = getword(sptr);

	if (!theword) {
		return (NULL);
	}
	theop = (op *) calloc(1, sizeof(op));

	/* check for MAGIC */
	if (!strcmp(theword, "MAGIC")) {
		theop->type = CONSTANT_TYPE;
		theop->opno = 0;
		theop->data = MAGIC;
		return (theop);
	}
	/* check for a constant number */
	if (1 == sscanf(theword, "%lg %n", &number, &nchars)) {
		if (nchars != strlen(theword)) {
			fprintf(stderr, "getop: can't understand %s\n", theword);
			exit(-1);
		} else {
			theop->type = CONSTANT_TYPE;
			theop->opno = 0;
			theop->data = number;
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
			if (1 == sscanf(theword, "%d %n", &fitsnum, &nchars)) {
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



