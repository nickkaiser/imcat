#include <stdio.h>

main (int argc, char *argv[])
{
	char	string[1024];
	int	nchars, res;
	double 	number;

	strcpy(string, "123.456");
	res = sscanf(string, "%lg%n", &number, &nchars);

	fprintf(stdout, "res=%d\nstring=%s\nnumber=%lf\nnchars=%d\n",
		res, string, number, nchars);

	exit(0);
}

