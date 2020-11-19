/*
 * ximfront.h
 */

/* database for xfgetwidgetnumber(), xfgetnumber()*/
#define TITLE		1
#define PRINT		2
#define MAKEMASK	3 
#define TOPLEVEL	4
#define	PICWIDGET	5
#define	POSLABEL	6
#define	PIXVALLABEL	7
#define	SLIDER2		8
#define PLOT		9

/* names for xfgetvalue() */
#define FMAX		100
#define FMIN		101
#define GOODSIZE	102

/* declarations */
void	Syntax(XtAppContext app_con, char *call);
void	makeslider(int index, Widget fromVert, Widget fromHoriz);
void	xfgetargstring(char **argstrptr);
void	makewidgetdatabase(void);
int	xfgetvalue(int thename);
int	xfgetwidgetnumber(Widget button);
Widget	xfgetwidget(int thenumber);
void 	sliderChanged(Widget w, XtPointer closure, XtPointer call_data);
void	setslider(int index, float thumbpos, char *label, char *value);
void	setlabelstring(Widget w, char *thestring);
void	setlabelnumber(Widget w, float thenumber);







