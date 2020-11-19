/*
 * front.h
 */

/* database for xfgetwidgetnumber(), xfgetnumber()*/
#define TOPLEVEL	4
#define	PICWIDGET	5
#define	FLABEL	6
#define	TLABEL	7

/* names for xfgetvalue() */
#define FMAX		100
#define FMIN		101
#define GOODSIZE	102

/* declarations */
void	Syntax(XtAppContext app_con, char *call);
void	xfgetargstring(char **argstrptr);
void	makewidgetdatabase(void);
float	xfgetvalue(int thename);
int	xfgetwidgetnumber(Widget button);
Widget	xfgetwidget(int thenumber);
void	setlabelstring(Widget w, char *thestring);
void	setlabelnumber(Widget w, float thenumber);







