#include <stdio.h>
#include <stdlib.h>
#include "extravars.h"


main(int argc, char * argv[])
{
	char 	*iplist, *opdefs;
	int	thesize;

	if (argc != 2) {
		fprintf(stderr, "usage : extravars_test definitionstring\n");
		exit(-1);
	}
	parseextravars(argv[1], &thesize, &iplist, &opdefs);
	fprintf(stderr, "size=%d\niplist=%s\nopdefs=%s\n", thesize, iplist, opdefs);
	exit(0);
}

