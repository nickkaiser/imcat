#include <stdio.h>
#include "error.h"
#include "iostream.h"

main(int argc, char *argv[]) {
	iostream *ipstream, *opstream;
	char	line[1024];

	if (argc != 3) {
		error_exit("usage: tmp instream outstream\n");
	}
	ipstream = openiostream(argv[1], "r");
	opstream = openiostream(argv[2], "w");

	while(fgets(line, 1024, ipstream->f)) {
		fputs(line, opstream->f);
	}
}