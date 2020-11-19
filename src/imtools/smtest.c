/*
 * smtest.c
 */
 
#include <stdio.h>
#include <stdlib.h>
#include "../utils/smpopup.h"

main()
{
	smpopup(drawfn, cursorfn, "-g 256x256+100+100");
}

void	drawfn(void)
{
	float	thepointtype = 40.0;

	sm_limits(0.0, 10.0, 0.0, 10.0);
	sm_box(1, 2, 0, 0);
	sm_ptype(&thepointtype, 1);
	
}

void	cursorfn(float x, float y, int key)
{
	fprintf(stderr, "key = %d; x = %10.3e; y = %10.3e\n", key, x, y);
}

