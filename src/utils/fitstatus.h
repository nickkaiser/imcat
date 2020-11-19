/*
 * fitstatus.h
 *
 * defines the values for the status value
 * and declares functions for getting and setting it
 */

#define DET_NEG                 1
#define TRACE_NEG               2
#define TOO_MANY_ITERATIONS     4

int	getfitstatus(void);
void	setfitstatus(int thestatus);
