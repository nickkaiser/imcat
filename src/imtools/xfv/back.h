/*
 * back.h
 */


void	readimheader(int *n1, int *n2);
void	set_shades(void);
void	alloc_shades(void);
void	create_pixmap(void);
void	fill_pixmap(void);
void	makeXImage(void);
void	showpixmap(void);
int	myworkproc1(XtPointer client_data);	/* background task */
void    showpixmap(void);
void	setblocksize(int b);
void	setfrange(float ffmin, float ffmax);
void	setcolormapindex(int theindex);












