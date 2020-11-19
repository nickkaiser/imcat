#define usage "\n\n\n\
NAME\n\
	tablesort --- sort a table of numbers\n\
SYNOPSIS\n\
	tablesort	col\n\
DESCRIPTION\n\
		tablesort read lines of a table from stdin\n\
		lines beginning with \"#\" and empty lines are ignored\n\
		lines sorted in ascending order of colth column values.\n\
		negative column number for descending order\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n\n"

#include <stdio.h>
#include <math.h>

#define MAX_LINE_LENGTH 	4096
#define MAX_COLS		32
#define MAX_WORD_LENGTH		128

typedef struct thing {
	struct thing 	*left;
	struct thing 	*right;
	char	line[MAX_LINE_LENGTH];
	double	x;
	int	occupied;
} thing;

char	theline[MAX_LINE_LENGTH];

main(int argc, char *argv[])
{
	char	s[MAX_COLS][MAX_WORD_LENGTH];
	int	thecol, ncols, descending = 0;
	thing	*root;
	double	x;
	
	if (argc != 2) {
		fprintf(stderr, "%s", usage);
		exit(-1);
	}
	
	if (argv[1][1] == 'u') {
		fprintf(stderr, "%s", usage);
		exit(-1);
	} else {
		sscanf(argv[1], "%d", &thecol);
		if (thecol < 0) {
			thecol = -thecol;
			descending = 1;
		}
	}
	
	root = (thing *) calloc(1, sizeof(thing));
	while (fgets(theline, MAX_LINE_LENGTH, stdin)) {
		if (theline[0] == '#' || theline[0] == '\n') {
			continue;
		}
		ncols = sscanf(theline, 
"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s ", 
s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], 
s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31] 
		);
		if (thecol > ncols) {
			fprintf(stderr, "tablesort: too few columns in input\n");
			exit(0);
		}
		sscanf(s[thecol - 1], "%lf", &x);
		install(x, root);
	}
	if (descending)
		writedesc(root);
	else
		writeasc(root);
}


int	install(double	x, thing *thething)
{
	if (!(thething->occupied)) {
		thething->x = x;
		thething->occupied = 1;
		strcpy(thething->line, theline);
	} else {
		if (x > thething->x) {
			if (!(thething->right))
				thething->right = (thing *) calloc(1, sizeof(thing));
			install(x, thething->right);
		} else {
			if (!(thething->left))
				thething->left = (thing *) calloc(1, sizeof(thing));
			install(x, thething->left);
		}
	}
}


int	writeasc(thing *thething)
{
	if (thething->left)
		writeasc(thething->left);
	fprintf(stdout, "%s", thething->line);
	if (thething->right)
		writeasc(thething->right);
}



int	writedesc(thing *thething)
{
	if (thething->right)
		writedesc(thething->right);
	fprintf(stdout, "%s", thething->line);
	if (thething->left)
		writedesc(thething->left);
}

