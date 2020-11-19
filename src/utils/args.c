/*
 * args.c
 *
 * Nick Kaiser 1998
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "args.h"

static int	garg, gargc;
static char	**gargv, *gusage;

/* initialise globals garg, gargc, gragv, gusage */
int	argsinit(int argc, char *argv[], char *theusage)
{
	garg = 1;
	gargc = argc;
	gargv = argv;
	gusage = theusage;
}

/* check next arg exists and begins with '-' and return the rest of the string */
char	*getflag(void)
{
	char	errstring[64];

	if (garg < gargc) {
		if (strncmp(gargv[garg], "-", 1)) {
			sprintf(errstring, "expected '-' prefixed option flag, found '%s'", gargv[garg]);
			argserror(errstring, 0);
		} else {
			if (strlen(gargv[garg]) < 2) {
				argserror("zero length flag", 0);
			} else {
				return (gargv[garg++]+1);
			}
		}
	} else {
		return NULL;
	}
}

void	argserror(char *errorstring, int argshift)
{
	if (argshift == -1) {
		fprintf(stderr, "%s : error after argument %d '%s' : %s\n", gargv[0], garg + argshift, gargv[garg + argshift], errorstring);
	} else {
		fprintf(stderr, "%s : error parsing argument %d '%s' : %s\n", gargv[0], garg + argshift, gargv[garg + argshift], errorstring);
	}
	if (gusage) {
		fprintf(stderr, "%s", gusage);
	}
	exit(1);
}

/* get argument string */
char	*getargs(void)
{
	if (garg < gargc) {
		return (gargv[garg++]);
	} else {
		fprintf(stderr, "getargs: insufficient arguments\n");
		if (gusage) {
			fprintf(stderr, "%s", gusage);
		}
		exit(1);
	}
}

/* get integer argument */
int	getargi(void)
{
	int	i, pos;

	if (garg < gargc) {
		if (1 == sscanf(gargv[garg], "%d%n", &i, &pos)) {
			if (pos != strlen(gargv[garg])) {
				argserror("unexpected characters in integer argument", 0);
			} else {
				garg++;
				return (i);
			}
		} else {
			argserror("expected integer argument", 0);
		}
	} else {
		argserror("expected integer argument", -1);
	}
}

/* get float argument */
float	getargf(void)
{
	int	pos;
	float	f;

	if (garg < gargc && 1 == sscanf(gargv[garg], "%f%n", &f, &pos)) {
		if (pos != strlen(gargv[garg])) {
			argserror("unexpected characters in floating point argument", 0);
		} else {
			garg++;
			return (f);
		}
	} else {
		argserror("expected floating point argument", (garg < gargc ? 0 : -1));
	}
}


/* get double argument */
double	getargd(void)
{
	int	pos;
	double	f;

	if (garg < gargc && 1 == sscanf(gargv[garg], "%lf%n", &f, &pos)) {
		if (pos != strlen(gargv[garg])) {
			argserror("unexpected characters in floating point argument", 0);
		} else {
			garg++;
			return (f);
		}
	} else {
		argserror("expected floating point argument", (garg < gargc ? 0 : -1));
	}
}

/* test value of next argument */
int	nextargtype(void)
{
	int	i, pos;
	float	f;

	if (garg < gargc) {
		if (1 == sscanf(gargv[garg], "%d%n", &i, &pos)) {
			if (pos == strlen(gargv[garg])) {
				return(INUM_ARG);
			}
		}
		if (1 == sscanf(gargv[garg], "%f%n", &f, &pos)) {
			if (pos == strlen(gargv[garg])) {
				return(FNUM_ARG);
			}
		}
		if (gargv[garg][0] == '-') {
			if (strlen(gargv[garg]) == 1) {
				return(STDSTREAM_ARG);
			} else {
				return(FLAG_ARG);
			}
		} else {
			return(TEXT_ARG);
		}		
	} else {
		return(NO_ARG);
	}
}
