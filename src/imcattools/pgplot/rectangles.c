/*
 * rectangles.c
 */

#include <stdlib.h>
#include <math.h>
#include "rectangles.h"

#define MAX_RECTS 1000

typedef struct rectangle {
	float x1, x2, y1, y2;
	int	active;
} rectangle;

static rectangle *rectlis[MAX_RECTS];
static int	pos = 0, lastpos;


int	installrect(float x1, float x2, float y1, float y2)
{
	rectlis[pos] = (rectangle *) calloc(1, sizeof(rectangle));
	if (x2 > x1) {
		rectlis[pos]->x1 = x1;
		rectlis[pos]->x2 = x2;
	} else {
		rectlis[pos]->x1 = x2;
		rectlis[pos]->x2 = x1;
	}
	if (y2 > y1) {
		rectlis[pos]->y1 = y1;
		rectlis[pos]->y2 = y2;
	} else {
		rectlis[pos]->y1 = y2;
		rectlis[pos]->y2 = y1;
	}
	rectlis[pos]->active = 1;
	pos++;
	lastpos = pos;
}


int	deactivate(float x, float y, float *x1, float *x2, float *y1, float *y2)
{
	float	X1, X2, Y1, Y2;
	int	active;

	lastpos = pos;
	pos = 0;

	while (pos < lastpos) {
		X1 = rectlis[pos]->x1;
		X2 = rectlis[pos]->x2;
		Y1 = rectlis[pos]->y1;
		Y2 = rectlis[pos]->y2;
		active = rectlis[pos]->active;
		if (active && X1 < x && x < X2 && Y1 < y && y < Y2) {
			rectlis[pos]->active = 0;
			*x1 = X1;
			*x2 = X2;
			*y1 = Y1;
			*y2 = Y2;
			pos = lastpos;
			return (1);
		}
		pos++;
	}
	pos = lastpos;
	return (0);
}


void	setposfirst(void)
{
	pos = 0;
}


int	getrect(float *x1, float *x2, float *y1, float *y2)
{
	if (pos < lastpos) {
		while (!(rectlis[pos]->active)) {
			pos++;
			if (pos == lastpos)
				return(0);
		}
		*x1 = rectlis[pos]->x1;
		*x2 = rectlis[pos]->x2;
		*y1 = rectlis[pos]->y1;
		*y2 = rectlis[pos]->y2;
		pos++;
		return (1);
	} else {
		return (0);
	}
}

