/*
 * function to parse a extravars definition string and return the overall size
 * in doubles and strings to be used in the input and output pipe commands
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void parseextravars(char *definitions, int *thesize, char **iplist, char **opdefs)
{
	int 	deflen, tokensize, first;
	char	*thename, *thesizestring, *thetoken, *tmpstring;

	deflen = strlen(definitions);
	*iplist = calloc(3 * deflen, sizeof(char));
	*opdefs = calloc(3 * deflen, sizeof(char));
	tmpstring = calloc(128, sizeof(char));
	
	*thesize = 0;
	first = 1;
	thename = strtok(definitions, ":");
	while(thetoken = strtok(NULL, ":")) {
		if (!first) {
			thename = thetoken;
			thesizestring = strtok(NULL, ":");
		} else {
			thesizestring = thetoken;
		}
		sprintf(tmpstring, "%s ", thename);
		strcat(*iplist, tmpstring);
		sscanf(thesizestring, "%d", &tokensize);
		*thesize += tokensize;
		if (tokensize > 1) {
			sprintf(tmpstring, "-N '1 %d %s' ", tokensize, thename);
		} else {
			sprintf(tmpstring, "-n %s ", thename);
		}
		strcat(*opdefs, tmpstring);
		first = 0;
	}
}
