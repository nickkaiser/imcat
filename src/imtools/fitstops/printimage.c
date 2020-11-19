/*
 * printimage.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "psutils.h"

void	printimage(int Ncolors, int N1, int N2, unsigned char ***c, int width, int height) 
{
	int	psstringlen;
	char	*psstring;
	FILE	*opf;
	int	i, j, color;

	opf = get_print_opf();	
	psstringlen = 2 * Ncolors * N1 + 2;
	psstring = (char *) calloc(psstringlen, sizeof(char));
	fprintf(opf, "/picstr %d string def\n", N1 * Ncolors);
	ps("0 0 translate");
	
	fprintf(opf, "%d %d scale\n", (int) width, (int) height);
	switch (Ncolors) {
		case 1:
			sprintf(psstring, 
				"/drawImage {%d %d 8 [%d 0 0 %d 0 0] {currentfile picstr readhexstring pop} image} def", 
				N1, N2, N1, N2);
			break;
		case 3:
			sprintf(psstring, 
				"/drawImage {%d %d 8 [%d 0 0 %d 0 0] {currentfile picstr readhexstring pop} false 3 colorimage} def", 
				N1, N2, N1, N2);
			break;
		default:
			error_exit("print_im: I only know how to deal with one or three colors\n");
			break;
	}
	ps(psstring);
	ps("drawImage");
	for (i = 0; i < N2; i++) {
		for(j = 0; j < N1; j++) {
			for (color = 0; color < Ncolors; color++) {
				if (c[color][i][j] < 16) {
					fprintf(opf, "0%X", (int) c[color][i][j]);
				} else {
					fprintf(opf, "%X", (int) c[color][i][j]);
				}
			}
		}
		fprintf(opf, "\n");
	}
}