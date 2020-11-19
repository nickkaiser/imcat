#include <stdio.h>
#include <stdlib.h>

main(int argc, char *argv[])
{
	FILE	*ipf;
	char	line[1024];

	ipf = popen("imhead -g 32 3 32 32 100", "r");
	fgets(line, 1024, ipf);
	printf("%s\n", line);
	pclose(ipf);
	exit(0);
}

