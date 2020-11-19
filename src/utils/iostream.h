/*
 * iostream.h - utilities for opening files, pipes
 *
 * programs which use these routines to open/close input/output streams
 * will transparently recognise '-' as indicating stdin or stdout
 *
 * they also recognise the syntax
 *	'command |' and '| command' for piping from/to a command
 */

/* structure to define a stream */
typedef struct iostream {
	FILE 	*f;
	char	*mode;
	int	type;
} iostream;

/* definition of various stream types */
#define	STD_IOSTREAM_TYPE	0
#define FILE_IOSTREAM_TYPE	1
#define PIPE_IOSTREAM_TYPE	2


iostream *openiostream(char *iostreamstring, char *mode);
int	closeiostream(iostream *theiostream);
