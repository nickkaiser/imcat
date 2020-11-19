/*
 * args.h
 *
 * Nick Kaiser 1998
 */

#define	CURRENT_ARG 	0
#define PREVIOUS_ARG	-1

#define NO_ARG 		0
#define FLAG_ARG	1
#define INUM_ARG	2
#define FNUM_ARG	3
#define TEXT_ARG	4
#define STDSTREAM_ARG	5

/* initialise globals garg, gargc, gragv */
int	argsinit(int argc, char *argv[], char *theusage);

/* check next arg exists and begins with '-' and return the rest of the string */
char	*getflag(void);

/* print progname, string and exit with non-zero error status */
void	argserror(char *errorstring, int argshift);

/* get argument string */
char	*getargs(void);

/* get integer argument */
int	getargi(void);

/* get float argument */
float	getargf(void);

/* get double argument */
double	getargd(void);

/* test value of next argument */
int	nextargtype(void);
