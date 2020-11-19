/*
 * rshloop.c 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "utils/error.h"
#include "utils/args.h"
#include "imlib/fits.h"

#define usage "\n\n\
NAME\n\
	rshloop - execute shell commands in parallel.\n\
\n\
SYNOPSIS\n\
	rshloop [-u] [-d] [-v] [-q] [-s slavelist] [-f ifmt] [-i initcom] ni commandformatstring\n\
\n\
DESCRIPTION\n\
	'rshloop' executes shell commands in parallel across a network.\n\
	It is currently functionally equivalent, and has identical arguments,\n\
	to 'pvmloop'.\n\
\n\
	'rshloop' first reads a table of slaves from 'rshloopslaves.lst' (or\n\
	from 'slavelist' with the -s option). The first line of\n\
	this table should contain ns, the number of slaves to be used,\n\
	followed by ns lines containing two strings: nodetag and nodename.\n\
	For example:\n\
\n\
		3\n\
		01	node01\n\
		02	node02\n\
		03	node03\n\
\n\
	You may specify any node, including the local node, multiple times.  This will\n\
	result in multiple processes running on each physical node.\n\
\n\
	'rshloop' then creates a set of temporary FIFO pipes /tmp/pid.s.fifo, one per slave.\n\
\n\
	'rshloop' then initiates a sequence of nb batches of nc <= ns commands until\n\
	all of the ni commands have been exhausted.\n\
	The form of the command executed is\n\
\n\
		rsh -n remotenode 'runcom command s' > /tmp/pid.s.fifo &\n\
\n\
	where 'command' is generated from the 'commandformatstring' as\n\
	described below and s is the the slave number.\n\
\n\
	If the environment variable 'RSHLOOP_RSH' is defined then it will\n\
	be used in place of 'rsh -n', so you can use e.g. ssh if you like.\n\
\n\
	The command 'runcom' executes the command given as its first argument\n\
	and saves the stdout and stderr, if not already redirected elsewhere to\n\
	temporary files. Once the command has terminated, it sends a message to\n\
	its stdout containing any stderr and stdout from the command and then\n\
	cleans up after itself.\n\
\n\
	Once a batch of commands has been initiated, rshloop opens all of the\n\
	FIFOs for reading (otherwise the commands will block).  It then loops\n\
	over the nc processes and collects the output.\n\
\n\
	The commandformatstring has a syntax similar to a printf format\n\
	string. On each iteration it is processed and each occurence\n\
	of %%i is replaced by the iteration number i = 0...ni-1,\n\
	%%I is replaced by a fixed length textual representation of i\n\
	%%Jn.m is replaced by a representation of i in base m with n digits\n\
	(n and m must both be single digit integers).\n\
	%%n is replaced by the node number, and %%N is replaced by the\n\
	nodetag, and %%%% is replaced by %%.  The result\n\
	may be a compound command of subprocesses linked by the\n\
	pipe symbol '|', and the output of the command may be\n\
	redirected into a disk file on the slave (in which case the master process\n\
	will receive, and generate, no standard output).\n\
\n\
EXAMPLES\n\
	These examples assume the rshloopslaves.lst file above, and that each\n\
	slave has a disk named /dnn, where nn is the nodetag.\n\
\n\
	To check on the status of the slaves:\n\
		rshloop 3 w\n\
\n\
	To clear a scratch directory on each slave:\n\
		rshloop 3 'rm /d%%N/tmp/*'\n\
\n\
	To generate a set of 1000 Monte Carlo simulations with some command 'monty'\n\
	which takes as an argument a seed (given here by the iteration number):\n\
		rshloop 1000 \"monty -seed %%i > /d%%N/tmp/monty%%I.dat\"\n\
	This would cause the following commands to be executed:\n\
		monty -seed 0   > /d01/tmp/monty000.dat   (on node01)\n\
		monty -seed 1   > /d02/tmp/monty001.dat   (on node02)\n\
		monty -seed 2   > /d03/tmp/monty002.dat   (on node03)\n\
		monty -seed 3   > /d01/tmp/monty003.dat   (on node01)\n\
		......\n\
		monty -seed 998 > /d01/tmp/monty998.dat   (on node03)\n\
		monty -seed 999 > /d01/tmp/monty999.dat   (on node01)\n\
\n\
OPTIONS\n\
	With -u flag we output this man page and exit.\n\
	With -d flag we just output the series of commands that would otherwise be executed\n\
	and exit. This is highly recommended with commands that delete files etc.\n\
	Use -v flag to invoke verbose mode, or -q to run quietly.\n\
	Use -f ifmt to specify a format string for the iteration number,\n\
	otherwise we use '%%.nd' specification where n is just large\n\
	enough to hold the number ni-1.\n\
	With '-i inticom' option the shell command string 'initcom' is executed\n\
	on the master node using system(initcom) before any other output is collected.\n\
	This allows you to generate a header, to which the output\n\
	of the slave processes can be prepended.  For example,\n\
	rshloop -q -i 'lc -C -n x < /dev/null' 10 'makerandcat 1000 -seed 2 -dim 1 | lc -o'\n\
	generates a catalogue containing 10 x 1000 random numbers.  The effect of the 'lc -o'\n\
	call is to chop the headers off the slave output.\n\
	Similarly you can use 'imhead -g ....' to generate a fits header and pipe\n\
	the output of the slaves through 'imhead -d'.\n\
\n\
BUGS\n\
	'rshloop' does not die if the remote shell commands fail.\n\
\n\
	'rshloop' does not do any load balancing.\n\
\n\
FILES\n\
	rshloopslaves.lst\n\
\n\
SEE ALSO\n\
	pvmloop, pvmserver, topvm.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n"

void	makecommand(char *srccom, char *dstcom, int i, char *ifmt, int s, char *nodetag);

#define MAX_SLAVES		10000
#define MAX_COMMAND_LEN		100000
#define MAX_NODENAME_LEN	1024
#define MAX_LINE_LEN		100000

#define OP_NULL		0
#define	OP_STDERR	1
#define OP_STDOUT	2
#define OP_FAILED	3
#define OP_SUCCESS	4

main(int argc, char *argv[])
{
	int	s, ns, i, ni, b, nb, c, nc, debugmode, verbosity, ilen, pid, op_mode;
	char	*flag, *srccom, *thecom, *fullcom, *tmpcom, *slavelistname, line[MAX_LINE_LEN], *logname;
	char	*ifmt, *initcom, *nodetag[MAX_SLAVES], *nodename[MAX_SLAVES], pipename[64], rshcom[1024];
	FILE	*rshloopdbf, *thepipe[MAX_SLAVES], *comfilepipe;

	char	*command;

	/* defaults */
	slavelistname = "rshloopslaves.lst";
	debugmode = 	0;
	verbosity = 	1;
	ifmt = 		NULL;
	initcom = 	NULL;

	/* get the RSHLOOP_RSH string if defined */
	if (getenv("RSHLOOP_RSH")) {
		strcpy(rshcom, getenv("RSHLOOP_RSH"));
	} else {
		strcpy(rshcom, "rsh -n");
	}

	/* parse args */
	argsinit(argc, argv, usage);
	while (nextargtype() == FLAG_ARG) {
		flag = getflag();
		switch (flag[0]) {
			case 'd':
				debugmode = 1;
				break;
			case 'v':
				verbosity = 2;
				break;
			case 'q':
				verbosity = 0;
				break;
			case 's':
				slavelistname = getargs();
				break;
			case 'f':
				ifmt = getargs();
				break;
			case 'i':
				initcom = getargs();
				break;
			case 'u':
			default:
				error_exit(usage);
				break;
		}
	}
	ni = getargi();
	srccom = getargs();

	if (verbosity > 1) {
		fprintf(stderr, "# rshloop: using remote command '%s'\n", rshcom);
	}

	if (!ifmt) {
		/* generate format string */
		ilen = (int) ceil(log10((double) ni));
		ifmt = (char *) calloc(8, sizeof(char));
		sprintf(ifmt, "%%.%dd", ilen);
	}

	/* allocate the command strings */
	thecom  = (char *) calloc(MAX_COMMAND_LEN, sizeof(char));
	fullcom = (char *) calloc(MAX_COMMAND_LEN, sizeof(char));
	tmpcom = (char *) calloc(MAX_COMMAND_LEN, sizeof(char));

	/* read slavelist file */
	rshloopdbf = fopen(slavelistname, "r");
	if (rshloopdbf) {
		if (1 != fscanf(rshloopdbf, "%d", &ns)) {
			error_exit("rshloop: slavelist format error\n");
		}
		for (s = 0; s < ns; s++) {
			nodetag[s] = (char *) calloc(64, sizeof(char));
			nodename[s] = (char *) calloc(64, sizeof(char));
			if (2 != fscanf(rshloopdbf, "%s %s", nodetag[s], nodename[s])) {
				error_exit("rshloop: slavelist format error\n");
			}
		}
	} else {
		error_exit("rshloop: failed to open slavelist file\n");
	}

	if (verbosity > 1) {
		/* echo a list of slave numbers etc to stderr */
		fprintf(stderr, "# rshloop: nslaves = %d:\n", ns);
		for (s = 0; s < ns; s++) {
			fprintf(stderr, "#\t%d\t%s\t%s\n", s, nodetag[s], nodename[s]);
		}
	}

	if (debugmode) {
		/* write the commands that would be executed to stdout */
		for (i = 0; i < ni; i++) {
			s = i % ns;
			makecommand(srccom, thecom, i, ifmt, s, nodetag[s]);
			/* sprintf the command to convert %% to % */
			sprintf(fullcom, thecom);
			fprintf(stdout, "%s: %s\n", nodename[s], fullcom);
		}
		exit(0);
	}

	if (initcom) {
		/* perform the initialisation command */
		if (system(initcom)) {
			fprintf(stderr, "rshloop: initialisation command ' %s ' failed\n", initcom);
			exit(-1);
		}
	}

	/* figure out how many batches we need to do */
	nb = (ni - ni % ns) / ns;
	if ((nb * ns) < ni) {
		nb++;
	}
	if (verbosity > 1) {
		fprintf(stderr, "# rshloop : ni=%d ns=%d nb=%d\n", ni, ns, nb);
	}

	/* get the process id */
	pid = (int) getpid();

	/* get the username */
	logname = getenv("LOGNAME");

	/* clean up any leftover tmp files */
	sprintf(tmpcom, "mkdir -p /tmp/rshloop_%s ; /bin/rm -f /tmp/rshloop_%s/runcom.%d.*", logname, logname, pid);
	if (verbosity > 1) {
		fprintf(stderr, "# rshloop : executing '%s'\n", tmpcom);
	}
	if (system(tmpcom)) {
		fprintf(stderr, "rshloop : command '%s' failed\n", tmpcom);
		exit(-1);
	}

	/* create the FIFOs  and also create the remote tmp directories */
	for (s = 0; s < ns; s++) {
		sprintf(tmpcom, "mkdir -p /tmp/rshloop_%s ; mknod /tmp/rshloop_%s/runcom.%d.%d.fifo p", logname, logname, pid, s);
		if (verbosity > 1) {
			fprintf(stderr, "# rshloop : executing '%s'\n", tmpcom);
		}
		if (system(tmpcom)) {
			fprintf(stderr, "rshloop : command '%s' failed\n", tmpcom);
			exit(-1);
		}
		sprintf(tmpcom, "%s %s 'mkdir -p /tmp/rshloop_%s'", rshcom, nodename[s], logname);
                if (verbosity > 1) {
                        fprintf(stderr, "# rshloop : executing '%s'\n", tmpcom);
                }
                if (system(tmpcom)) {
                        fprintf(stderr, "rshloop : command '%s' failed\n", tmpcom);
                        exit(-1);
                }
	}

	/* loop over batches */
	for (b = 0; b < nb; b++) {
		nc = ni - b * ns;
		nc = (nc > ns ? ns : nc);
		if (verbosity > 1) {
			fprintf(stderr, "# rshloop : starting batch of %d processes\n", nc);
		}
		/* start the processes */
		for (c = 0; c < nc; c++) {
			i = b * ns + c;
			/* generate the command */
			makecommand(srccom, thecom, i, ifmt, c, nodetag[c]);
			/* output message to stderr */
			fprintf(stderr, "# slave %d : %s\n", c, thecom);
			/* create the command to write the command file */
			if (strcmp(rshcom, "rsh -n")) {
				sprintf(tmpcom, "%s %s 'cat > /tmp/rshloop_%s/runcom.%d.%d.com'", rshcom, nodename[c], logname, pid, c);
			} else {
				sprintf(tmpcom, "rsh %s 'cat > /tmp/rshloop_%s/runcom.%d.%d.com'", nodename[c], logname, pid, c);
			}
			comfilepipe = popen(tmpcom, "w");
			if (!comfilepipe) {
				fprintf(stderr, "rshloop : failed to open pipe to write command file\n");
				exit(1);
			}
			fprintf(comfilepipe, thecom);
			pclose(comfilepipe);
			if (verbosity > 1) {
				sprintf(fullcom, "%s %s 'runcom /tmp/rshloop_%s/runcom.%d.%d.com -v' > /tmp/rshloop_%s/runcom.%d.%d.fifo &", rshcom, nodename[c], logname, pid, c, logname, pid, c);
			} else {
				sprintf(fullcom, "%s %s 'runcom /tmp/rshloop_%s/runcom.%d.%d.com' > /tmp/rshloop_%s/runcom.%d.%d.fifo &", rshcom, nodename[c], logname, pid, c, logname, pid, c);
			}
			if (verbosity > 1) {
				fprintf(stderr, "# rshloop : executing command '%s'\n", fullcom);
			}
			if (system(fullcom)) {
				fprintf(stderr, "# rshloop : command '%s' failed\n", fullcom);
				exit(1);
			} else {
				if (verbosity > 1) {
					fprintf(stderr, "# rshloop : command '%s' successfully spawned\n", fullcom);
				}
			}
		}
		/* open the FIFOs - otherwise processes will block */
		for (c = 0; c < nc; c++) {
			sprintf(pipename, "/tmp/rshloop_%s/runcom.%d.%d.fifo", logname, pid, c);
			if (verbosity > 1) {
				fprintf(stderr, "# rshloop : opening pipe '%s'\n", pipename);
			}
			thepipe[c] = fopen(pipename, "r");
			if (thepipe){
				if (verbosity > 1) {
					fprintf(stderr, "# rshloop : successfully opened pipe '%s'\n", pipename);
				}
			} else {
				fprintf(stderr, "rshloop : failed to open pipe '%s'\n", pipename);
				exit(1);
			}
		}
		/* read the output and close the FIFOs */
		for (c = 0; c < nc; c++) {
			op_mode = OP_NULL;
			while(fgets(line, MAX_LINE_LEN, thepipe[c])) {
				if (!(strncmp(line, "# runcom : stdout", 17))) {
					op_mode = OP_NULL;
				}
				if (!(strncmp(line, "# runcom : all done", 19))) {
					op_mode = OP_NULL;
				}
				if (!(strncmp(line, "# runcom : failed", 17))) {
					op_mode = OP_FAILED;
				}
				if (!(strncmp(line, "# runcom : success", 18))) {
					op_mode = OP_SUCCESS;
				}
				switch (op_mode) {
					case OP_STDERR:
						fprintf(stderr, "rshloop : stderr from slave %d : %s", c, line);
						break;
					case OP_STDOUT:
						fprintf(stdout, "%s", line);
						break;
					case OP_FAILED:
						fprintf(stderr, "rshloop : slave %d failed : %s", c, line);
						break;
					case OP_SUCCESS:
						if (verbosity > 1) {
							fprintf(stderr, "rshloop : slave %d success : %s", c, line);
						}
						break;
				}
				if (!(strncmp(line, "# runcom : stderr", 17))) {
					op_mode = OP_STDERR;
				}
				if (!(strncmp(line, "# runcom : stdout", 17))) {
					op_mode = OP_STDOUT;
				}
			}
			fclose(thepipe[c]);
		}
	} 

		

	/* remove tmp files */
	sprintf(tmpcom, "/bin/rm -f /tmp/rshloop_%s/runcom.%d.*", logname, pid);
	if (verbosity > 1) {
		fprintf(stderr, "# rshloop : executing '%s'\n", tmpcom);
	}
		if (system(tmpcom)) {
			fprintf(stderr, "rshloop : command '%s' failed\n", tmpcom);
			exit(-1);
		}

	/* all done */
	if (verbosity > 1) {
		fprintf(stderr, "# rshloop: all done...\n");
	}
	exit(0);

}

void	makecommand(char *srccom, char *dstcom, int i, char *ifmt, int s, char *nodetag)
{
	char	tmp[4], digit[11] = "0123456789";
	int 	n, m, j, ii;

	while (1) {
		if (*srccom == '%') {
			srccom++;
			switch (*srccom) {
				case 'i':
					sprintf(dstcom, "%d", i);
					break;
				case 'n':
					sprintf(dstcom, "%d", s);
					break;
				case 'I':
					sprintf(dstcom, ifmt, i);
					break;
				case 'N':
					sprintf(dstcom, "%s", nodetag);
					break;
				case 'J':
					if (strlen(srccom) < 4) {
						error_exit("rshloop: command string format error\n");
					}
					strncpy(tmp, srccom + 1, 3);
					tmp[3] = '\0';
					if (2 != sscanf(tmp, "%d.%d", &n, &m)) {
						error_exit("rshloop: command string format error\n");
					}
					ii = i;
					for (j = 0; j < n; j++) {
						dstcom[n - j - 1] = digit[ii % m];
						ii = (ii - (ii % m)) / m;
					}
					dstcom[n] = '\0';
					srccom += 3;
					break;
				case '%':
					sprintf(dstcom, "%s", "%%");
					break;
				default:
					error_exit("rshloop: command string format error\n");
			}
			dstcom += strlen(dstcom);
			srccom++;
		} else {
			*dstcom = *srccom;
			if (*srccom == '\0') {
				return;
			}
			dstcom++;
			srccom++;
		}
	}
}



