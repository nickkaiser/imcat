#define usage "\n\n\n\
NAME\n\
	tablecumulate --- cumulate a column of numbers in a table\n\
\n\
SYNOPSIS\n\
	tablecumulate	col\n\
DESCRIPTION\n\
		tablecumulate read lines of a table containing lines\n\
			X_1 X_2 X_3 .....\n\
		from stdin\n\
		lines beginning with \"#\" and empty lines are ignored\n\
		cumulates the col'th column entry\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"

#include <stdio.h>
#include <math.h>

#define MAX_LINE_LENGTH 	4096
#define MAX_COLS		32
#define MAX_WORD_LENGTH		128

main(int argc, char *argv[])
{
	double 	x, xlast = 0.0;
	char	line[MAX_LINE_LENGTH], s[MAX_COLS][MAX_WORD_LENGTH];
	int		col, thecol, ncols;
	
	if (argc != 2) {
		fprintf(stderr, "%s", usage);
		exit(-1);
	}
	
	if (argv[1][1] == 'u') {
		fprintf(stderr, "%s", usage);
		exit(-1);
	}
	
	sscanf(argv[1], "%d", &thecol);
	
	fprintf(stdout, "# %s %s\n", argv[0], argv[1]);
	while (fgets(line, MAX_LINE_LENGTH, stdin)) {
		if (line[0] == '#' || line[0] == '\n') {
			continue;
		}
		ncols = sscanf(line, 
"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s ", 
s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], 
s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31] 
		);
		if (thecol > ncols) {
			fprintf(stderr, "tablecumulate: too few columns in input\n");
			exit(0);
		}
		sscanf(s[thecol - 1], "%lf", &x);
		xlast += x;
		sprintf(s[thecol - 1], "%g", xlast);
		for (col = 0; col < ncols; col++) {
			fprintf(stdout, "%s\t", s[col]);
		}
		fprintf(stdout, "\n");
	}
}

