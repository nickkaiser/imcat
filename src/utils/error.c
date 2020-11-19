/**
 ** miscellaneous service routines
 **/
 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "error.h"

void	error_exit(char *message)
{
	fprintf(stderr, message);
	exit(1);
}


