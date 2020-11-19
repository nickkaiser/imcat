#define usage "\n\n\n\
NAME\n\
	tableop --- apply operation to table of numbers\n\
SYNOPSIS\n\
	tableop	col op\n\
DESCRIPTION\n\
		tableop read lines of a table containing lines\n\
			X_1 X_2 X_3 .....\n\
		from stdin\n\
		lines beginning with \"#\" and empty lines are ignored\n\
		changes the col'th column entry\n\
		op can be one of	exp, ln, dex, log\n\
		e.g. tableop 2 ln will replace 2nd column its natural log\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"

#include <stdio.h>
#include <math.h>

#define MAX_LINE_LENGTH 	4096
#define MAX_COLS		32
#define MAX_WORD_LENGTH		128
#define	EXP	0
#define	LN	1
#define DEX	3
#define	LOG	4

main(int argc, char *argv[])
{
	double 	x, operand;
	char	line[MAX_LINE_LENGTH], s[MAX_COLS][MAX_WORD_LENGTH];
	int	col, thecol, ncols, op;
	
	if (argc != 3) {
		fprintf(stderr, "%s", usage);
		exit(-1);
	} else {
		sscanf(argv[1], "%d", &thecol);
		switch (argv[2][1]) {
			case 'x':
				op = EXP;
				break;
			case 'n':
				op = LN;
				break;
			case 'e':
				op = DEX;
				break;
			case 'o':
				op = LOG;
				break;
			default:
				fprintf(stderr, "%s", usage);
				exit(-1);
				break;
		}
	}
	
	fprintf(stdout, "# %s %s %s\n", argv[0], argv[1], argv[2]);
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
			fprintf(stderr, "tableop: too few columns in input\n");
			exit(0);
		}
		sscanf(s[thecol - 1], "%lf", &x);
		switch (op) {
			case EXP:
				x = exp(x);
				break;
			case LN:
				x = (x > 0 ? log(x) : 0);
				break;
			case DEX:
				x = pow(10.0, x);
				break;
			case LOG:
				x = (x > 0 ? log10(x) : 0);
				break;
		}
		sprintf(s[thecol - 1], "%g", x);
		for (col = 0; col < ncols; col++) {
			fprintf(stdout, "%s\t", s[col]);
		}
		fprintf(stdout, "\n");
	}
}

