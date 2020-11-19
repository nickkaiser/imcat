#define usage "\n\n\n\
NAME\n\
	tablesum --- sum a table of numbers\n\
\n\
SYNOPSIS\n\
	tablesum c1 p1 [c2 p2 ....] \n\
DESCRIPTION\n\
		tablesum read lines of a table containing lines\n\
			X_1 X_2 X_3 .....\n\
		from stdin\n\
		lines beginning with \"#\" and empty lines are ignored\n\
		returns sum X_c1^p1\n\
		if additional c,p pairs are given the it returns sum of\n\
		(X_c1^p1) * ( X_c2^p2) *.....\n\
		max line length 4096 characters, max number of cols = 32.\n\
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
	double 	x, sum, increment;
	char	line[MAX_LINE_LENGTH], s[MAX_COLS][MAX_WORD_LENGTH];
	int	col, ncols, c[MAX_COLS];
	double	p[MAX_COLS];
	
	if (argc < 3 || 2 * (argc / 2) == argc) {
		fprintf(stderr, "%s", usage);
		exit(-1);
	} else {
		ncols = (argc - 1) / 2;
		for (col = 0; col < ncols; col++) {
			sscanf(argv[2 * col+1], "%d", c+col);
			sscanf(argv[2 * col+2], "%lf", p+col);
			c[col]--;
		}
	}
	
	sum = 0.0;
	while (fgets(line, MAX_LINE_LENGTH, stdin)) {
		if (line[0] == '#' || line[0] == '\n') {
			continue;
		}
		sscanf(line, 
"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s ", 
s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], 
s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31] 
		);
		increment = 1.0;
		for (col = 0; col < ncols; col++) {
			if (1 != sscanf(s[c[col]], "%lf", &x)) {
				fprintf(stderr, "tablesum: failed to convert input\n");
				exit(-1);
			}
			increment *= pow(x, p[col]);
		}
		sum += increment;
	}
	fprintf(stdout, "%lf\n", sum);
}

