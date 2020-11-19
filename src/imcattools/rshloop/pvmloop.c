/*
 * pvmloop.c 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "pvm3.h"
#include "utils/error.h"
#include "utils/args.h"
#include "imlib/fits.h"
#include "pvmpack.h"
#include "pvmspawn.h"
#include "pvmtags.h"

#define usage "\n\n\
NAME\n\
	pvmloop - execute shell commands across a parallel virtual machine running PVM.\n\
\n\
SYNOPSIS\n\
	pvmloop [-u] [-d] [-v] [-q] [-s slavelist] [-f ifmt] [-i initcom] niter commandformatstring\n\
\n\
DESCRIPTION\n\
	'pvmloop' executes shell commands across a parallel virtual machine.\n\
\n\
	'pvmloop' first attempts to enroll in pvm, assumed to be running,\n\
	and issues an error message and exits if it fails. See\n\
	http://www.epm.ornl.gov/pvm/pvm_home.html for more information about PVM.\n\
\n\
	'pvmloop' then reads a table of slaves from 'pvmslaves.lst' (or\n\
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
	'pvmloop' then spawns a set of 'pvmserver' proceses, one per slave,\n\
	which are then send a sequence of shell commands of the form:\n\
\n\
		( ( command ) | topvm ...) 2>&1 topvm ...\n\
\n\
	where 'command' is generated from the 'commandformatstring' as\n\
	described below.  The purpose of the 'topvm' processes is to\n\
	send the stdout and stderr of 'command' back to the master\n\
	process as pvm messages where they are decoded and merged into\n\
	the master process stdout and stderr streams.\n\
	By default, each complete command is echoed to stderr, but you\n\
	can switch this off with the -q option.\n\
\n\
	The commandformatstring has a syntax similar to a printf format\n\
	string. On each iteration it is processed and each occurence\n\
	of %%i is replaced by the iteration number i = 0...ni-1,\n\
	%%I is replaced by a fixed length textual representation of i\n\
	%%n is replaced by the node number, and %%N is replaced by the\n\
	nodetag, and %%%% is replaced by %%.  The result\n\
	may be a compound command of subprocesses linked by the\n\
	pipe symbol '|', and the output of the command may be\n\
	redirected into a disk file on the slave (in which case the master process\n\
	will receive, and generate, no standard output).\n\
	The command thus generated is then sent to a 'pvmserver'\n\
	to be executed.\n\
\n\
EXAMPLES\n\
	These examples assume the pvmslaves.lst file above, and that each\n\
	slave has a disk named /dnn, where nn is the nodetag.\n\
\n\
	To check on the status of the slaves:\n\
		pvmloop 3 w\n\
\n\
	To clear a scratch directory on each slave:\n\
		pvmloop 3 'rm /d%%N/tmp/*'\n\
\n\
	To generate a set of 1000 Monte Carlo simulations with some command 'monty'\n\
	which takes as an argument a seed (given here by the iteration number):\n\
		pvmloop 1000 \"monty -seed %%i > /d%%N/tmp/monty%%I.dat\"\n\
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
	pvmloop -q -i 'lc -C -n x < /dev/null' 10 'makerandcat 1000 -seed 2 -dim 1 | lc -o'\n\
	generates a catalogue containing 10 x 1000 random numbers.  The effect of the 'lc -o'\n\
	call is to chop the headers off the slave output.\n\
	Similarly you can use 'imhead -g ....' to generate a fits header and pipe\n\
	the output of the slaves through 'imhead -d'.\n\
\n\
BUGS\n\
	'pvmloop' does not die if the remote shell commands fail.\n\
\n\
	'pvmloop' does not do any load balancing.\n\
\n\
FILES\n\
	pvmslaves.lst\n\
\n\
SEE ALSO\n\
	pvmserver, topvm, apply.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n"

#define MAX_ARGS        128
#define MAX_ARG_LEN     1024
	
void	makefullcommand(char *srccom, char *dstcom, int mtid, int opmsgtag, int errmsgtag);
void	makecommand(char *srccom, char *dstcom, int i, char *ifmt, int s, char *nodetag);

main(int argc, char *argv[])
{
	int	theslave, thetid, mytid, dsttid, s, ns, tid[MAX_SLAVES];
	char	*flag, *srccom, *thecom, *fullcom, *slavelistname, *pvmargv[MAX_ARGS];
	char	*nodetag[MAX_SLAVES], *nodename[MAX_SLAVES], *buff, *ifmt, *initcom;
	FILE	*pvmdbf;
	int	commsgtag, outputmsgtag[MAX_SLAVES], errormsgtag[MAX_SLAVES];
	int	i, ni, debugmode, pvmarg, pvmargc, buffsize, ilen, status;
	int     ihost, nhosts, narch, verbosity, ntask, count;
	struct 	pvmhostinfo     *hostp = NULL;
	struct	pvmtaskinfo	*taskp = NULL;

	/* defaults */
	slavelistname = "pvmslaves.lst";
	debugmode = 	0;
	verbosity = 	1;
	ifmt = 		NULL;
	initcom = 	NULL;

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

	if (!ifmt) {
		/* generate format string */
		ilen = (int) ceil(log10((double) ni));
		ifmt = (char *) calloc(8, sizeof(char));
		sprintf(ifmt, "%%.%dd", ilen);
	}

	/* allocate the command strings */
	thecom  = (char *) calloc(100000, sizeof(char));
	fullcom = (char *) calloc(100000, sizeof(char));

	/* enroll in pvm and get our tid */
	mytid = pvm_mytid();
	if (mytid < 0) {
		error_exit("pvmloop: failed to enroll in pvm -- try starting pvm.\n");
	}
	pvm_config(&nhosts, &narch, &hostp);
	if (verbosity > 1) {
		fprintf(stderr, "# pvmloop: master tid = %x\n", mytid);
	}

	/* read slavelist file */
	pvmdbf = fopen(slavelistname, "r");
	if (pvmdbf) {
		if (1 != fscanf(pvmdbf, "%d", &ns)) {
			error_exit("pvmloop: slavelist format error\n");
		}
		for (s = 0; s < ns; s++) {
			nodetag[s] = (char *) calloc(64, sizeof(char));
			nodename[s] = (char *) calloc(64, sizeof(char));
			if (2 != fscanf(pvmdbf, "%s %s", nodetag[s], nodename[s])) {
				error_exit("pvmloop: slavelist format error\n");
			}
		}
	} else {
		error_exit("pvmloop: failed to open slavelist file\n");
	}

	if (verbosity > 1) {
		fprintf(stderr, "# pvmloop: nslaves = %d:\n", ns);
		for (s = 0; s < ns; s++) {
			fprintf(stderr, "#\t%d\t%s\t%s\n", s, nodetag[s], nodename[s]);
		}
	}

	/* make the message tags */
	for(s = 0; s < ns; s++) {
		outputmsgtag[s] = getuniquemsgtag(s);
		errormsgtag[s] = getuniquemsgtag(TAG_MULTIPLIER / 2 + s);
	}
	commsgtag = getuniquemsgtag(TAG_MULTIPLIER - 1);

	if (debugmode) {
		for (i = 0; i < ni; i++) {
			s = i % ns;
			makecommand(srccom, thecom, i, ifmt, s, nodetag[s]);
			fprintf(stdout, "%s: %s\n", nodename[s], thecom);
		}
		pvm_exit();
		exit(0);
	}

	if (initcom) {
		/* perform the initialisation command */
		if (system(initcom)) {
			error_exit("pvmloop: initialisation command ' %s ' failed\n");
		}
	}

	/* spawn one pvmserver processe per node */
	pvmargc = 2;
	for (pvmarg = 0; pvmarg < pvmargc; pvmarg++) {
		pvmargv[pvmarg] = (char *) calloc(16, sizeof(char));
	}
	pvmargv[pvmargc] = NULL;
	sprintf(pvmargv[0], "%x", mytid);
	sprintf(pvmargv[1], "%d", commsgtag);
	for (s = 0; s < ns; s++) {
		if (1 != pvm_spawn("pvmserver", pvmargv, PvmTaskHost, nodename[s], 1, tid + s)) {
			pvm_exit();
			error_exit("pvmloop: failed to spawn pvmserver process\n");
		} else {
			if (verbosity > 1) {
				fprintf(stderr, "# pvmloop: spawned 'pvmserver ");
				for (pvmarg = 0; pvmarg < pvmargc; pvmarg++) {
					fprintf(stderr, "%s ", pvmargv[pvmarg]);
				}
				fprintf(stderr, "' with tid %x on node %s\n", tid[s], nodename[s]);
			}
		}
	}


	/* send a set of commands, one per slave */
	for (i = 0; i < ni + ns; i++) {
		s = i % ns;
		if (i >= ns) {
			/* receive the output */
			pvm_recv(-1, outputmsgtag[s]);
			if (verbosity > 1) {
				fprintf(stderr, "# pvmloop: output %d received from node%s\n", i - ns, nodetag[s]);
			}
			buff = NULL;
			while (buffsize = unpacksegment(&buff)) {
				if (buffsize > 0) {
					fwrite(buff, sizeof(char), buffsize, stdout); 
					free(buff);
				}
			}
			buff = NULL;
			/* receive the error */
			pvm_recv(-1, errormsgtag[s]);
			buff = NULL;
			while (buffsize = unpacksegment(&buff)) {
				if (buffsize > 0) {
					fprintf(stderr, "# pvmloop: stderr from node%s: ", nodetag[s]);
					fwrite(buff, sizeof(char), buffsize, stderr); 
					free(buff);
				}
			}
			buff = NULL;
		}
		if (i < ni) {
			/* send a command */
			makecommand(srccom, thecom, i, ifmt, s, nodetag[s]);
			makefullcommand(thecom, fullcom, mytid, outputmsgtag[s], errormsgtag[s]);
			switch (verbosity) {
				case 0:
					break;
				case 1:
					fprintf(stderr, "# pvmloop: %s: command: %s\n", nodename[s], thecom);
					break;
				default:
					fprintf(stderr, "# pvmloop: %s: command: %s\n", nodename[s], fullcom);
					break;				
			}
			/* pack the command */
			pvm_initsend(PvmDataDefault);
			packstring("command:");
			packstring(fullcom);
			/* and send it */
			pvm_send(tid[s], commsgtag);
		}
	}

	/* all done */
	if (verbosity > 1) {
		fprintf(stderr, "# pvmloop: killing slaves...\n");
	}
	killslaves(ns, tid, commsgtag);
	if (verbosity > 1) {
		fprintf(stderr, "# pvmloop: exiting pvm...\n");
	}
	pvm_exit();
	if (verbosity > 1) {
		fprintf(stderr, "# pvmloop: all done...\n");
	}
	exit(0);

}

void	makecommand(char *srccom, char *dstcom, int i, char *ifmt, int s, char *nodetag)
{
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
				case '%':
					sprintf(dstcom, "%s", "%%");
					break;
				default:
					error_exit("pvmloop: command string format error\n");
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

void	makefullcommand(char *srccom, char *dstcom, int mtid, int opmsgtag, int errmsgtag)
{
	sprintf(dstcom, "(( ");
	dstcom += strlen(dstcom);
	sprintf(dstcom, srccom);
	dstcom += strlen(dstcom);
	sprintf(dstcom, " | topvm %x %d ) 2>&1 | topvm %x %d )", mtid, opmsgtag, mtid, errmsgtag);
	return;
}


