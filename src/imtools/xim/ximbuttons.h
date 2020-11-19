/*
 * ximbuttons.h
 */

void	Zoom(Widget w, XtPointer client_data, XtPointer call_data);
void	UnZoom(Widget w, XtPointer client_data, XtPointer call_data);
void	makepopup(Widget button, XtPointer client_data, XtPointer call_data);
void	Cancel(Widget button, XtPointer client_data, XtPointer call_data);
void	OkPrint(Widget button, XtPointer client_data, XtPointer call_data);
void	OkMakemask(Widget button, XtPointer client_data, XtPointer call_data);
void 	Makemask(void);
void	Print(void);
void	printclick(int i, int j);
void	Title(Widget w, XtPointer client_data, XtPointer call_data);
void	picclick(Widget w, XtPointer client_data, XtPointer call_data);
void	openmaskfile(char *maskfilename);
void	closemaskfile(void);
void	addmaskpoint(int i, int j);
void	Contour(Widget button, XtPointer client_data, XtPointer call_data);
void	Profile(Widget button, XtPointer client_data, XtPointer call_data);
void	Custom(Widget button, XtPointer client_data, XtPointer call_data);
void	plotclick(int i, int j);
void	doplot(int icentre, int jcentre, int iside, int jside, int level, int plottype);






