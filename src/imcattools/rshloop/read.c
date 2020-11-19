#include <stdio.h>

#define usage "usage: read fifo\n"

#define BUFFSIZE	1024

main(int argc, char *argv[])
{
	char	*fifoname, ipbuff[BUFFSIZE];
	FILE	*ipf;

	if (argc != 2) {
		fprintf(stderr, usage);
		exit(1);
	} else {
		fifoname = argv[1];
	}

	fprintf(stderr, "# read: opening fifo %s\n", fifoname);
	ipf = fopen(fifoname, "r");
	if (ipf) {
		fprintf(stderr, "# read: fifo %s opened successfully\n", fifoname);
	} else {
		fprintf(stderr, "# read: failed to open fifo %s\n", fifoname);
		exit(1);
	}
	while(fgets(ipbuff, BUFFSIZE, ipf)) {
		fprintf(stdout, ipbuff);
	}

	fprintf(stderr, "# read: exhausted fifo %s\n", fifoname);
	exit(0);
}
