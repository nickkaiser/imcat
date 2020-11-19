#include <stdio.h>
#include <unistd.h>

#define usage "\nNAME\n\
	sleeper - sleep and then return stdout and stderr test messages\n\
\n\
SYNOPSIS\n\
	sleeper processid delay dostdout dostderr\n\
\n\
DESCRIPTION\n\
	sleeper sleeps for delay seconds.  Then, if dostdout is\n\
	non zero, it sends to stdout the message\n\
		# stdout from procid processid\n\
	If dostderr is non-zero it sends to stderr, the string\n\
		# stderr from procid processid\n\
	If dostderr is greater than 1 it exits with failure\n\
	status, otherwise it returns with success status.\n\
\n\
SEE ALSO\n\
	rshloop runcom\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"


main(int argc, char *argv[])
{
	int	procid, dostdout, dostderr, delay;

	if (argc != 5) {
		fprintf(stderr, usage);
		exit(1);
	} else {
		sscanf(argv[1], "%d", &procid);
		sscanf(argv[2], "%d", &delay);
		sscanf(argv[3], "%d", &dostdout);
		sscanf(argv[4], "%d", &dostderr);
	}

	sleep(delay);
	if (dostdout) {
		fprintf(stdout, "# stdout from procid %d\n", procid);
	}
	if (dostderr) {
		fprintf(stderr, "# stderr from procid %d\n", procid);
	}
	if (dostderr > 1) {
		exit(-1);
	} else {
		exit(0);
	}
}
