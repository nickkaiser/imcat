/*
 * psutils.h
 */

void	ps(char *string);
void	psDrawChar(char aChar);
void	print_caption(char *caption);
void	set_print_opf(FILE *thefile);
FILE	*get_print_opf(void);
void	psline(double x1, double y1, double x2, double y2);
