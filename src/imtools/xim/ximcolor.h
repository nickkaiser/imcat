/*
 * ximcolor.h
 */

#define Q0_MAX	1
#define A_MAX	4

void	setcolorscheme(void);
void	set_shades(void);
int	color_index(int level, int i, int j);
float	newslidervalue(int index);
void	printall(int level, char *psfilename);
void	printselection(int level, int i1, int j1, int i2, int j2, 
		char *psfilename);
void	colorsetlevel(int level);
double	gray(int i, int j, int color);
void	setsliderval(int index, double theval);















