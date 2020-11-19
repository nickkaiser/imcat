#include <stdio.h>

main(int argc, char *argv[])
{
	char	command[1024], comfilename[64], line[1024], pipename[64];
	int	i, nproc;
	FILE	*comfile, *thepipe[2];

	nproc = 2;

	fprintf(stderr, "# testrshloop - testing runcom\n");
	/* start the processes */
	for (i = 0; i < nproc; i++) {
		sprintf(comfilename, "/tmp/tmp.d.com", i);
		comfile = fopen(comfilename, "w");
		/* generate a command to sleep for 2 seconds with procid = i and to return both stdin and stdout */
		fprintf(comfile, "sleeper %d 2 1 1", i);
		fclose(comfile);
		sprintf(command, "runcom %s", comfilename);
		fprintf(stderr, "# testrshloop : executing command '%s'\n", command);
		if (system(command)) {
			fprintf(stderr, "# testrshloop : command '%s' failed\n");
			exit(1);
		} else {
			fprintf(stderr, "# testrshloop : command '%s' successfully spawned\n");
		}
	}

exit(0);

	for (i = 0; i < nproc; i++) {
		sprintf(pipename, "stdout%d", i);
		fprintf(stderr, "# testrshloop : opening pipe '%s'\n", pipename);
		thepipe[i] = fopen(pipename, "r");
		if (thepipe[i]){
			fprintf(stderr, "# testrshloop : successfully opened pipe '%s'\n", pipename);
		} else {
			fprintf(stderr, "# testrshloop : failed to open pipe '%s'\n", pipename);
			exit(1);
		}
	}

	for (i = 0; i < nproc; i++) {
		while(fgets(line, 1024, thepipe[i])) {
			fprintf(stderr, line);
		}
		fclose(thepipe[i]);
	}
	exit(0);
}
