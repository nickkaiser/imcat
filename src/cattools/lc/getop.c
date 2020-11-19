/*
 * getop.c
 *
 */

#include <stdio.h>
#include "../../catlib/cat.h"
#include "getop.h"
#include "operators.h"
#include "error.h"

static	cathead	*sourcecathead;
static	cathead	*destcathead;


item	*getop(char **sptr)
{
	item	*theitem, *sourceitem;
	char	*theword, *wordbase;
	double	number;
	int	nchars;
	
	wordbase = theword = getword(sptr);
	if (!theword) {
		return (NULL);
	}
	/* first we check for a string */
	if (theword[0] == '{') {
		if (theword[strlen(theword) - 1] != '}') {
			error_exit("getop: unmatched '{'\n");
		}
		theword++;
		theword[strlen(theword) - 1] = '\0';
		theitem = newitem("VAR", TEXT_TYPE, 1, 1);
		allocitemcontents(theitem, &(theitem->addr), 0);
		copystring(&(((char **)(theitem->addr))[0]), theword);
		return (theitem);
	}
	/* now we check for a number */
	/* used to have "%lg %n" here, but it broke on the mac */
	if (1 == sscanf(theword, "%lg%n", &number, &nchars)) {
		if (nchars != strlen(theword)) {
			fprintf(stderr, "getop: can't understand %s\n", theword);
			exit(-1);
		} else {
			theitem = newitem("VAR", NUM_TYPE, 1, 1);
			allocitemcontents(theitem, &(theitem->addr), 0);
			((double *)(theitem->addr))[0] = number;
			return (theitem);
		}
	}
	/* now check for variables, operators */
	switch (theword[0]) {
		case '%':			/* on object item */
			theword++;
			if (!strlen(theword)) {
				error_exit("getop: expression format error\n");
			}
			sourceitem = getobjectitem(theword, sourcecathead);
			theitem = copyitem(sourceitem);
			theitem->addr = sourceitem->addr;
			break;
		case '^':			/* a header item */
			theword++;
			if (!strlen(theword)) {
				error_exit("getop: expression format error\n");
			}
			if (!(sourceitem = getheaderitem(theword, sourcecathead))) {
				error_exit("getop: missing necessary header item\n");
			}
			theitem = copyitem(sourceitem);
			theitem->addr = sourceitem->addr;
			break;
		default:
			theitem = newopitem(theword);
			break;		
	}
	free (wordbase);
	return(theitem);
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



void	setsourcecathead(cathead *thecathead)
{
	sourcecathead = thecathead;
}


void	setdestcathead(cathead *thecathead)
{
	destcathead = thecathead;
}
