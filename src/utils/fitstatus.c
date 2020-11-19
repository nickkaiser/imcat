/*
 * fitstatus.c
 *
 * defines the values for the status value
 * and provides functions for getting and setting it
 */

#define DET_NEG                 1
#define TRACE_NEG               2
#define TOO_MANY_ITERATIONS     4

static int 	fitstatus;

int	getfitstatus(void)
{
	return(fitstatus);
}

void	setfitstatus(int thestatus)
{
	fitstatus = thestatus;
}
