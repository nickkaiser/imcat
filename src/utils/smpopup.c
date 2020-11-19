/*
 * smpopup.c
 */
 
/*
 * function to
 *	0) send a message to stderr describing smpopup key interpretation
 *	1) pop up a sm window (with specified geometry if geomstring is not NULL)
 *	2) call the drawfn() to draw window contents
 *	3) put up a cursor and handle key/mouse inputs as follows:
 *		if key = r we recall drawfn() and then call cursorfn()
 *		if key = x we quit the popup and then call cursorfn()
 *		otherwise we just call cursorfn()
 */

#include <stdio.h>
#include <stdlib.h>
#include "smpopup.h"

void	smpopup(void (*drawfn)(void), void (*cursorfn)(float x, float y, int key), char *geomstring)
{
	float	x, y;
	int		key;
	char	devicestring[128];
	
	/* print the key interpretation message */
	fprintf(stderr, "# sm-popup: key = 'r' to redraw, key = 'x' to exit popup\n");
	strcpy(devicestring, "X11 ");
	strcat(devicestring, geomstring);
	sm_device(devicestring);
	sm_graphics();
	sm_erase();
	drawfn();
	sm_gflush();
	while(1) {
		sm_curs(&x, &y, &key);
		switch (key) {
			case 'r':
				drawfn();
				sm_gflush();
				break;
			case 'x':
				return;
				break;
			default:
				break;
		}
		cursorfn(x, y, key);
	}
}


