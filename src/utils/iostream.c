/*
 * iostream.h - utilities for opening files, pipes
 */

#include <stdio.h>
#include "iostream.h"
#include "error.h"

int hasbadchars(char *string);

iostream *openiostream(char *iostreamstring, char *mode)
{
	iostream	*thestream;
	int		i;

	/* allocate the iostream structure */
	thestream = (iostream *) calloc(1, sizeof(iostream));
	if (!thestream) {
		error_exit("openiostream: failed to allocate the iostream structure\n");
	}
	
	/* analyse the mode */
	if (strlen(mode) != 1 || (mode[0] != 'r' && mode[0] != 'w')) {
		error_exit("openiostream: illegal mode argument\n");
	}
	thestream->mode = (char *) calloc(2, sizeof(char));
	strcpy(thestream->mode, mode);
		
	/* analyse the iostreamstring */
	if (strlen(iostreamstring) < 1) {
		error_exit("openiostream: empty iostreamstring argument\n");
	}
	if (iostreamstring[0] == '-') {
		if (strlen(iostreamstring) > 1) {
			fprintf(stderr, "openiostream: warning filename %s begins with hyphen\n", iostreamstring);
		} else {
			thestream->type = STD_IOSTREAM_TYPE;
			switch (mode[0]) {
				case 'r':
					thestream->f = stdin;
					return (thestream);
				case 'w':
					thestream->f = stdout;
					return (thestream);
			}
		}
	}
	switch (mode[0]) {
		case 'r':
			if (iostreamstring[strlen(iostreamstring) - 1] == '|') {
				thestream->type = PIPE_IOSTREAM_TYPE;
				iostreamstring[strlen(iostreamstring) - 1] = '\0';
				while (iostreamstring[strlen(iostreamstring) - 1] == ' ') {
					iostreamstring[strlen(iostreamstring) - 1] = '\0';
				}
				thestream->f = popen(iostreamstring, "r");
			} else {
				thestream->type = FILE_IOSTREAM_TYPE;
				if (hasbadchars(iostreamstring)) {
					fprintf(stderr, "openiostream: filename %s contains bad characters\n", iostreamstring);
					exit(1);
				}
				thestream->f = fopen(iostreamstring, "r");
			}
			break;
		case 'w':
			if (iostreamstring[0] == '|') {
				thestream->type = PIPE_IOSTREAM_TYPE;
				iostreamstring++;
				while (iostreamstring[0] == ' ') {
					iostreamstring++;
				}
				thestream->f = popen(iostreamstring, "w");
			} else {
				thestream->type = FILE_IOSTREAM_TYPE;
				if (hasbadchars(iostreamstring)) {
					fprintf(stderr, "openiostream: filename '%s' contains bad characters\n", iostreamstring);
					exit(1);
				}
				thestream->f = fopen(iostreamstring, "w");
			}
			break;
	}
	if (thestream->f == NULL) {
		error_exit("openiostream: failed to open stream\n");
	}		
	return (thestream);
}

int	closeiostream(iostream *theiostream)
{
	switch(theiostream->type) {
		case STD_IOSTREAM_TYPE:
			break;
		case FILE_IOSTREAM_TYPE:
			fclose(theiostream->f);
			break;
		case PIPE_IOSTREAM_TYPE:
			pclose(theiostream->f);
			break;
		default:
			error_exit("closeiostream: illegal type\n");
	}
	free(theiostream);
}


int	hasbadchars(char *string)
{
	int 	i;

	for (i = 0; i < strlen(string); i++) {
		if (string[i] == ' ' || string[i] == '|') {
			return(1);
		}
	}
	return(0);
}
