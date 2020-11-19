/*
 * ximback.h
 */

#define	MAX_LEVELS	10
#define MAX_COLORS 	3


/* names for xbgetvalue(), xbsetvalue() */
#define LEVEL		0
#define NLEVELS		1
#define NCOLORS		2
#define LASTUSERI	3
#define LASTUSERJ	4

int	xbgetvalue(int thename);
void	xbsetvalue(int thename, int thevalue);
void	alloc_shades(void);
Pixmap	create_pixmap(int level);
void	fill_pixmap(void);
void	makepreview(void);
void	readdataheader(void);
void	allocatedata(void);
void	readdata(void);
short	**alloc_f(int N1, int N2);
void	makescrunchedviews(void);
void	makeXImage(int N1, int N2);
void	GetPixelValue(int pixi, int pixj, int height, int width,
			int *useri, int *userj, short *pixval);
Pixmap	createzoompixmap(void);
void	fillzoompixmap(void);
void	destroyzoompixmap(void);
void	showpixmap(void);
int	myworkproc1(XtPointer client_data);	/* background task */
int	myworkproc2(XtPointer client_data);	/* background task */
int	myworkproc3(XtPointer client_data);	/* background task */











