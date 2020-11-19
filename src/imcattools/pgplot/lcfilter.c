/*
 * lcfilter.c
 */

#include <stdio.h>
#include "lcfilter.h"

#define	LCFILTERSTRINGLENGTH 4096

static 	char 	lcfilter[LCFILTERSTRINGLENGTH], xname[64], yname[64];
static	int	firstclause, mode;

void	startfilter(int coordsargcount, char *coordsarg[], int editmode)
{
	switch (coordsargcount) {
		case 1:
			sprintf(xname, "%%%s[0]", coordsarg[0]);
			sprintf(yname, "%%%s[1]", coordsarg[0]);
			break;
		case 2:
			sprintf(xname, "%%%s", coordsarg[0]);
			sprintf(yname, "%%%s", coordsarg[1]);
			break;
		default:
			fprintf(stderr, "plotcat: invalid coordsargcount\n");
			exit(-1);
			break;
	}
	mode = editmode;
	firstclause = 1;
	sprintf(lcfilter, "lc -b -i '");
}


int	addfiltercondition(float x1, float x2, float y1, float y2)
{
	char	clause[512];

	sprintf(clause, "%s %g > %s %g < and %s %g > and %s %g < and ",
		xname, x1, xname, x2, yname, y1, yname, y2);
	if (!firstclause)
		strcat(clause, "or ");
	if (strlen(lcfilter) + strlen(clause) > LCFILTERSTRINGLENGTH - 30) {
		fprintf(stderr, "plotcat: filter string too long!\n");
		exit(-1);
	}
	firstclause = 0;
	strcat(lcfilter, clause);		
}


void	dofilter(char *tempfilename)
{
	FILE	*filterpipe;

	if (mode == REJECT_MODE) {
	strcat(lcfilter, "!");
	}
	strcat(lcfilter, "' < ");
	strcat(lcfilter, tempfilename);
	
	if (!(filterpipe = popen(lcfilter, "w"))) {
		fprintf(stderr, "plotcat: failed to open filter pipe!\n");
		exit(-1);
	}
	pclose(filterpipe);	
}
