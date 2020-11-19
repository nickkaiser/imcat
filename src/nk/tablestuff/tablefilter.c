#define usage "\n\n\n\
NAME\n\
	tablefilter --- filter a table of numbers\n\
\n\
SYNOPSIS\n\
	tablefilter col minval maxval\n\
DESCRIPTION\n\
		tablefilter read lines from a table from stdin\n\
		lines beginning with \"#\" and empty lines are passed as comments\n\
		as are lines in which value of the col'th entry\n\
		has (x >= minval && x <= maxval).\n\
		Entries spearated by tab(s) or space(s).  Ist column is col = 1\n\
		max line length 4096 characters\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"

#include <stdio.h>
#include <math.h>

#define MAX_LINE_LENGTH 	4096

main(int argc, char *argv[])
{
	double 	x, xmin, xmax;
	char	line[MAX_LINE_LENGTH];
	int	pos, lastwaswhite, col, thecol;
	
	if (argc != 4) {
		fprintf(stderr, "%s", usage);
		exit(-1);
	} else {
		sscanf(argv[1], "%d", &thecol);
		sscanf(argv[2], "%lf", &xmin);
		sscanf(argv[3], "%lf", &xmax);
	}
	
	fprintf(stdout, "# %s %s %s %s\n", argv[0], argv[1], argv[2], argv[3]);
	while (fgets(line, MAX_LINE_LENGTH, stdin)) {
		if (line[0] == '#' || line[0] == '\n') {
			fprintf(stdout, "%s", line);
			continue;
		}
		/* now we step over thecol - 1 columns */
		col = 0;
		lastwaswhite = 1;
		pos = 0;
		while (col < thecol) {
			if (line[pos] == ' ' || line[pos] == '\t') {
				lastwaswhite = 1;
				pos++;
				continue;
			} else {
				if (lastwaswhite)
					col++;
				lastwaswhite = 0;
				pos++;
			}
		}
		pos--;
		if (1 != sscanf(line + pos, "%lf", &x)) {
			fprintf(stderr, "tablefilter: failed to read an entry\n");
			exit(-1);
		} else {
			if (x >= xmin && x <= xmax)
				fprintf(stdout, "%s", line);
		}
	}
}

