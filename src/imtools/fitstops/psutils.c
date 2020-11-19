/*
 * psutils.c
 */

#include <stdio.h>
#include "psutils.h"


static char	psstring[10000];
static FILE	*opf;

void	ps(char *string)
{
	fprintf(opf, "%s\n", string);
}


void	psDrawChar(char aChar)
{
	switch (aChar) {
		case '(':
			sprintf(psstring, "(\\50) show");
			break;
		case ')':
			sprintf(psstring, "(\\51) show");
			break;
		default:
			sprintf(psstring, "(%c) show", aChar);
			break;
	}
	ps(psstring);
}



void	print_caption(char *caption) {
	float	size = 10, space = 13, h = 13, v = -28;
	int		c, len;
	len = strlen(caption);
	
	ps("/Times-Roman findfont   10.000 scalefont setfont");
	sprintf(psstring, "%f %f moveto", h, v);
	ps(psstring);
	for (c = 0; c < len; c++)
		switch (caption[c]) {
			case '\n':
				v -= space;
				sprintf(psstring, "%f %f moveto", h, v);
				ps(psstring);
				break;
			default:
				psDrawChar(caption[c]);
		}
}


void	set_print_opf(FILE *thefile)
{
	opf = thefile;
}

FILE	*get_print_opf(void)
{
	return (opf);
}


void	psline(double x1, double y1, double x2, double y2)
{
	fprintf(opf, "n %.2lf %.2lf m %.2lf %.2lf l\n", x1, y1, x2, y2); 
}
