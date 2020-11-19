#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>

#define usage "\nNAME\n\
        runcom - run a command for rshloop\n\
\n\
SYNOPSIS\n\
        runcom comfile [-v]\n\
\n\
DESCRIPTION\n\
        runcom reads a command from comfile and then executes it,\n\
	capturing the stderr and stdout in temporary files\n\
	/tmp/rshloop_username/runcom.PID.stderr and /tmp/rshloop_username/runcom.PID.stderr.\n\
	It then sends those to stdout.\n\
\n\
SEE ALSO\n\
        rshloop\n\
AUTHOR\n\
        Nick Kaiser --- kaiser@hawaii.edu\n\n"


main(int argc, char *argv[])
{
	char	*comfilename, thecommand[10000], fullcommand[10000], tmpcommand[1024], *logname;
	pid_t	pid;
	FILE	*comfile;
	int	verbose;

	verbose = 0;

	switch (argc) {
		case 3:
			verbose = 1;
		case 2:
			comfilename = argv[1];
			break;
		default:
			fprintf(stderr, usage);
			exit(1);
	}

	if (!strcmp(comfilename, "-u")) {
		fprintf(stderr, usage);
		exit(1);
	}

	if (!(comfile = fopen(comfilename, "r"))) {
		fprintf(stdout, "# runcom : failed to open %s\n", comfilename);
		exit(1);
	}
	fgets(thecommand, 10000, comfile);
	fclose(comfile);
	sprintf(tmpcommand, "/bin/rm -f %s", comfilename);
	system(tmpcommand);	

	pid = getpid();
	logname = getenv("LOGNAME");

	sprintf(fullcommand, "( ( %s ) > /tmp/rshloop_%s/runcom.%d.stdout ) 2> /tmp/rshloop_%s/runcom.%d.stderr", thecommand, logname, (int) pid, logname, (int) pid);

	if (system(fullcommand)) {
		fprintf(stdout, "# runcom : failed to execute command : %s\n", fullcommand);
	} else {
		if (verbose) {
			fprintf(stdout, "# runcom : successfully executed command : %s\n", fullcommand);
		}
	}
	fflush(stdout);
	fprintf(stdout, "# runcom : stderr:\n");
	fflush(stdout);
	sprintf(tmpcommand, "/bin/cat /tmp/rshloop_%s/runcom.%d.stderr", logname, (int) pid);
	system(tmpcommand);
	fflush(stdout);
	fprintf(stdout, "# runcom : stdout:\n");
	sprintf(tmpcommand, "/bin/cat /tmp/rshloop_%s/runcom.%d.stdout", logname, (int) pid);
	fflush(stdout);
	system(tmpcommand);	
	fflush(stdout);
	fprintf(stdout, "# runcom : all done\n");
	
	sprintf(tmpcommand, "/bin/rm -f /tmp/rshloop_%s/runcom.%d.stdout /tmp/rshloop_%s/runcom.%d.stderr", logname, (int) pid, logname, (int) pid);
	system(tmpcommand);
	exit(0);
}
