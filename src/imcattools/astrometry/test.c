#include <stdio.h>
#include <math.h>

main (int argc, char *argv[])
{
	char	format[64];
	int	i;

	
	fprintf(stdout, "             +\n");
	sprintf(format, "%s", argv[1]);
	for (i = -10; i < 10; i++) {
		fprintf(stdout, format, M_PI * pow(10.0, i));
		fprintf(stdout, "+\n");
	}
}
