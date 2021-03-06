#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#define usage "\nNAME\n\
        piepsend - send dummy data to a pipe\n\
\n\
SYNOPSIS\n\
        pipesend nblocks blocksize [-f outputfile]\n\
\n\
DESCRIPTION\n\
        pipesend sends blocks of dummy data to stdout.\n\
	It reports to stderr the number of bytes sent and the\n\
	elapsed time.\n\
\n\
	With the -f option it writes to a file with fwwrite and\n\
	then performs the fsync command.\n\
\n\
	It may be useful for testing network or disk performance.\n\
\n\
SEE ALSO\n\
	piperecv\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\n"

main (int argc, char *argv[]) {
	int 	nblocks, blocksize, i, j;
	long	nbytes;
	char	*buff[2], *filename;
	double	starttime, endtime, elapsedtime, Mbps;
	struct timeval	tv;
	struct timezone	tz;
	FILE	*opf;

	/* defaults */
	filename = NULL;
	opf = stdout;

	if (argc < 3) {
		fprintf(stderr, usage);
		exit(-1);
	}
	sscanf(argv[1], "%d", &nblocks);
	sscanf(argv[2], "%d", &blocksize);
	nbytes = nblocks * blocksize;
	if (argc == 5) {
		filename = argv[4];
		opf = fopen(filename, "w");
	}

	buff[0] = (char *) calloc(blocksize, sizeof(char));
	buff[1] = (char *) calloc(blocksize, sizeof(char));
	for (j = 0; j < blocksize; j++) {
		buff[0][j] = 'a';
		buff[1][j] = 'b';
	}

	gettimeofday(&tv, &tz);
	starttime = tv.tv_sec + 1.e-6 * tv.tv_usec;

	fprintf(opf, "%d %d\n", nblocks, blocksize);
	for (i = 0; i < nblocks; i++) {
		fwrite(buff[i % 2], sizeof(char), blocksize, opf);
	}

	if (filename) {
		fsync(fileno(opf));
	}

	gettimeofday(&tv, &tz);
	endtime = tv.tv_sec + 1.e-6 * tv.tv_usec;
	elapsedtime = endtime - starttime;
	Mbps = 8.e-6 * nbytes / elapsedtime;

	fprintf(stderr, "# pipesend : %ld bytes transferred in %14.8lf seconds: rate = %14.8lf Mbps\n", 
		nbytes, elapsedtime, Mbps);
	
	free(buff[0]);
	free(buff[1]);
	exit(0);
}
