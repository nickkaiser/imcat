#include <stdio.h>

#define usage "usage: maketestcat nobj ncols\n\
SUPERSEDED BY makerandcat and makegridcat\n"

double	drand48();

main(int argc, char *argv[]){
	int	nlines, ncols, line, col;

	if (argc != 3) {
		fprintf(stderr, usage);
		exit(-1);
	}
	sscanf(argv[1], "%d", &nlines);
	sscanf(argv[2], "%d", &ncols);

	for (line = 0; line < nlines; line++) {
		for (col = 0; col < ncols; col++) {
			printf("%13.8lg ", drand48());
		}
		printf("\n");
	}
}
