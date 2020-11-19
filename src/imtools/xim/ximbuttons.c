/*
 * ximbuttons.c
 */
 
#include 	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<limits.h>

#include <X11/Xlib.h>
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Shell.h>
#include <X11/cursorfont.h>
#include <X11/Xaw/Dialog.h>	

#include	"ximbuttons.h"
#include	"ximfront.h"
#include	"ximback.h"
#include	"../../imlib/fits.h"

static	int		makingmask = 0, printing = 0, gotfirstprintclick;
static	int		plotting = 0, gotfirstplotclick, plottype;
static	int		firstmaskpoint;
static	int		hitok, hitcancel;
static	char		filename[128], customcommand[256];
static	FILE		*maskf;

extern	int		N1[], N2[];
extern	short	**f[MAX_FITS_DIM][MAX_LEVELS];

#define      MAX(x,y) (((x) > (y)) ? (x) : (y)) 
#define      MIN(x,y) (((x) < (y)) ? (x) : (y)) 
#define      ABS(x) ((x < 0) ? (-x) : (x)) 

#define	CONTOURPLOT	0
#define PROFILEPLOT	1
#define CUSTOMPLOT	2

#define MAGIC SHORT_MAGIC

void	Zoom(Widget w, XtPointer client_data, XtPointer call_data)
{
	int	level;

	level = xbgetvalue(LEVEL);
	if (level == 0)
		createzoompixmap();
	level--;
	xbsetvalue(LEVEL, level);
      	if (level >= 0) {
		showpixmap();
	} else {
		fillzoompixmap();
		showpixmap();
	}
	
}

void	UnZoom(Widget w, XtPointer client_data, XtPointer call_data)
{
	int	level, nlevels;

	level = xbgetvalue(LEVEL);
	nlevels = xbgetvalue(NLEVELS);

	if (level < nlevels - 1) {
		level++;
		level = (level < 0 ? 0 : level);
		xbsetvalue(LEVEL, level);
 		if (level == 0)
			destroyzoompixmap();
		showpixmap();
	}
}


void	makepopup(Widget button, XtPointer client_data, XtPointer call_data)
{
	Arg			arg[5];
	Position	x, y;
	Dimension	height, width;
	Cardinal	n;
	Widget		popup, dialog;
	int		buttonnumber;
	char		*stringlist[] = {
"xim: X image viewer by Nick Kaiser\0",
"to create a postscript file:\n\
1) hit ok button\n\
2) select rectangle to print or\n\
   hit print button again to print\n\
   the entire view\0",
"to create a mask file:\n\
1) hit ok button\n\
2) click pairs of points to define rectangles\n\
3) hit makemask button again when you're done\0",
"to plot contours or profile of subimage:\n\
1) hit contour or profile button\n\
2) 2 clicks to define centre and edge of square subimage\0",
NULL};

	buttonnumber = xfgetwidgetnumber(button);

	if ((buttonnumber == PRINT) && printing) {		/* second hit on print button */
		Print();
		return;
	}
	if ((buttonnumber == MAKEMASK) && makingmask) {		/* second hit on makemask button */
		Makemask();
		return;
	}

	/* we want the credits dialog to come up over the title button */
	n = 0;
	XtSetArg(arg[0], XtNheight, &height); n++;
	XtSetArg(arg[1], XtNwidth, &width); n++;
	XtGetValues(button, arg, n);
	XtTranslateCoords(button, (Position) (width / 2), (Position) (height / 2),
		      &x, &y);
	n = 0;
	XtSetArg(arg[n], XtNx, x);				n++;
	XtSetArg(arg[n], XtNy, y);				n++;
	popup = XtCreatePopupShell("prompt", transientShellWidgetClass, button,
			       arg, n);
    	switch (buttonnumber) {
		case TITLE:
			dialog = XtVaCreateManagedWidget("dialog", dialogWidgetClass, popup,
				XtNlabel,	stringlist[0],
			NULL);
    			XawDialogAddButton(dialog, "ok", Cancel, (XtPointer) dialog);
			break;
		case PRINT:
			dialog = XtVaCreateManagedWidget("dialog", dialogWidgetClass, popup,
				XtNlabel,	stringlist[1],
				XtNvalue,	"xxxx.ps",
				NULL);
    				XawDialogAddButton(dialog, "ok", OkPrint, (XtPointer) dialog);
  				XawDialogAddButton(dialog, "cancel", Cancel, (XtPointer) dialog);
			break;
		case MAKEMASK:
			dialog = XtVaCreateManagedWidget("dialog", dialogWidgetClass, popup,
				XtNlabel,	stringlist[2],
				XtNvalue,	"xxxx.msk",
				NULL);
    				XawDialogAddButton(dialog, "ok", OkMakemask, (XtPointer) dialog);
  				XawDialogAddButton(dialog, "cancel", Cancel, (XtPointer) dialog);
			break;
		case PLOT:
			dialog = XtVaCreateManagedWidget("dialog", dialogWidgetClass, popup,
				XtNlabel,	stringlist[3],
				XtNvalue,	" ",
				NULL);
    				XawDialogAddButton(dialog, "contour", Contour, (XtPointer) dialog);
    				XawDialogAddButton(dialog, "profile", Profile, (XtPointer) dialog);
    				XawDialogAddButton(dialog, "custom string", Custom, (XtPointer) dialog);
  				XawDialogAddButton(dialog, "cancel", Cancel, (XtPointer) dialog);
			break;
		default:
			error_exit("makepopup: bad index\n");
	}
	XtPopup(popup, XtGrabNone);
}


void	Cancel(Widget button, XtPointer client_data, XtPointer call_data)
{
    Widget popup = XtParent( (Widget) client_data);

    hitcancel = 1;
    hitok = 0;
    XtDestroyWidget(popup);
}


void	OkPrint(Widget button, XtPointer client_data, XtPointer call_data)
{
	Widget popup = XtParent( (Widget) client_data);

	hitcancel = 0;
	hitok = 1;
	strcpy(filename, XawDialogGetValueString((Widget) client_data)); 
	XtDestroyWidget(popup);
	Print();
}


void	OkMakemask(Widget button, XtPointer client_data, XtPointer call_data)
{
	Widget popup = XtParent( (Widget) client_data);

	hitcancel = 0;
	hitok = 1;
	strcpy(filename, XawDialogGetValueString((Widget) client_data)); 
	XtDestroyWidget(popup);
	Makemask();
}



void	Contour(Widget button, XtPointer client_data, XtPointer call_data)
{
	Widget popup = XtParent( (Widget) client_data);

	plotting = 1;
	gotfirstplotclick = 0;
	plottype = CONTOURPLOT;
	XtDestroyWidget(popup);
}


void	Profile(Widget button, XtPointer client_data, XtPointer call_data)
{
	Widget popup = XtParent( (Widget) client_data);

	plotting = 1;
	gotfirstplotclick = 0;
	plottype = PROFILEPLOT;
	XtDestroyWidget(popup);
}



void	Custom(Widget button, XtPointer client_data, XtPointer call_data)
{
	Widget popup = XtParent( (Widget) client_data);

	plotting = 1;
	gotfirstplotclick = 0;
	plottype = CUSTOMPLOT;
	strcpy(customcommand, XawDialogGetValueString((Widget) client_data)); 
	XtDestroyWidget(popup);
}



void Makemask(void)
{
	if (!makingmask) {
		if (hitok) {
			makingmask = 1;
			openmaskfile(filename);
		}
	} else {
		makingmask = 0;
		closemaskfile();
	}
}



void	Print(void)
{
	int	level;

	level = xbgetvalue(LEVEL);

	if (!printing) {
		if (hitok) {
			printing = 1;
		}
	} else {
		printall((level < 0 ? 0 : level), filename);
		printing = gotfirstprintclick = 0;
	}
}


void	printclick(int i, int j)
{
	static	int	i1, i2, j1, j2, temp, level;
	
	level = xbgetvalue(LEVEL);

	if (!gotfirstprintclick) {
		i1 = i;
		j1 = j;
		gotfirstprintclick = 1;
	} else {
		i2 = i;
		j2 = j;
		printing = gotfirstprintclick = 0;
		if (i1 > i2) {
			temp = i1; i1 = i2; i2 = temp;
		}
		if (j1 > j2) {
			temp = j1; j1 = j2; j2 = temp;
		}
		printselection((level < 0 ? 0 : level), i1, j1, i2, j2, filename);
	}
}



void	plotclick(int i, int j)
{
	static	int	icentre, jcentre, iside, jside, level;
	
	if (gotfirstplotclick) {
		plotting = gotfirstplotclick = 0;
		iside = i;
		jside = j;
		level = xbgetvalue(LEVEL);
		level = (level < 0 ? 0 : level);
		doplot(icentre, jcentre, iside, jside, level, plottype);
	}  else {
		gotfirstplotclick = 1;
		icentre = i;
		jcentre = j;
	}
}


void	Title(Widget w, XtPointer client_data, XtPointer call_data)
{
	makepopup(w, client_data, call_data);
}




void	picclick(Widget w, XtPointer client_data, XtPointer call_data)
{
	Window		awindow, thewindow, root;
	Display		*thedisplay;
	int		i, j, pixi, pixj, useri, userj, ncolors;
	unsigned int	mask;
	Cardinal		n;
	Arg		arg[5];
	char	poslabelstring[20], pixvallabelstring[20];
	Dimension	width, height, borderwidth;
	short	pixval[MAX_COLORS];
	int		color;

	/* get the coordinates of pointer in window */
	thedisplay = XtDisplay(xfgetwidget(TOPLEVEL));
	thewindow = XtWindow(xfgetwidget(PICWIDGET));
	XQueryPointer(thedisplay, thewindow, &root,
		&awindow, &i, &j, &pixi, &pixj, &mask);
		
	/* get the border, height and width */
	n = 0;
    	XtSetArg(arg[n], XtNborderWidth, &borderwidth);	n++;
    	XtSetArg(arg[n], XtNheight, &height);	n++;
    	XtSetArg(arg[n], XtNwidth, &width);	n++;
	XtGetValues(xfgetwidget(PICWIDGET), arg, n);
		
	/* update the poslabel and pixvallabel */
	GetPixelValue(pixi, pixj, height, width, &useri, &userj, pixval);
	sprintf(poslabelstring, "i,j: %d %d\n",useri, userj);
	setlabelstring(xfgetwidget(POSLABEL), poslabelstring);
	ncolors = xbgetvalue(NCOLORS);
	switch(ncolors) {
		case 1:
			if (pixval[0] != MAGIC)
				sprintf(pixvallabelstring, "pixval: %d\n", (int) pixval[0]);
			else
				sprintf(pixvallabelstring, "pixval: MAGIC\n");
			break;
		case 2:
			sprintf(pixvallabelstring, "pixval: %d %d\n", (int) pixval[0], (int) pixval[1]);
			break;
		case 3:
			sprintf(pixvallabelstring, "pixval: %d %d %d\n", 
				(int) pixval[0], (int) pixval[1], (int) pixval[2]);
			break;
		default:
			break;
	}
	setlabelstring(xfgetwidget(PIXVALLABEL), pixvallabelstring);
	fprintf(stdout, "%4d %4d", useri, userj);
	for (color = 0; color < ncolors; color++)
		fprintf(stdout, " %6d", pixval[color]);
	fprintf(stdout, "\n");
	xbsetvalue(LASTUSERI, useri);
	xbsetvalue(LASTUSERJ, userj);
	if (makingmask)
		addmaskpoint(useri, userj);
	if (printing)
		printclick(useri, userj);
	if (plotting)
		plotclick(useri, userj);
}


void	openmaskfile(char *maskfilename)
{
	if (!(maskf = fopen(maskfilename, "w")))
		error_exit("openmaskfile: failed to open mask file for output\n");
	fprintf(stderr, "sending i,j coords to %s\n", maskfilename);
	fprintf(maskf, "#\tj1\ti1\tj2\ti2\n");
	firstmaskpoint = 1;
}



void	closemaskfile(void)
{
	fclose(maskf);
	fprintf(stderr, "closing mask file\n");
}




void	addmaskpoint(int i, int j)
{
	static	i1, i2, j1, j2, l, b, r, t;
	
	if (firstmaskpoint) {
		i1 = i;
		j1 = j;
		firstmaskpoint = 0;
	} else {
		i2 = i;
		j2 = j;
/*		if (i1 < i2) {
			b = i1;
			t = i2;
		} else {
			t = i1;
			b = i2;
		}
		if (j1 < j2) {
			l = j1;
			r = j2;
		} else {
			r = j1;
			j = j2;
		}
		fprintf(maskf, "%6d %6d %6d %6d\n", l, b, r, t);*/
		fprintf(maskf, "\t%6d\t%6d\t%6d\t%6d\n", j1, i1, j2, i2);
		firstmaskpoint = 1;
	}
}


void	doplot(int icentre, int jcentre, int iside, int jside, int level, int plottype)
{
	int		i, j, di, dj, io, jo, N;
	short	*ff;
	FILE	*fitspipe;
	fitsheader	*fits;
	fitscomment	*com;
	char	 astring[128], plotcommand[256];
	
	di = MAX(icentre, iside) - MIN(icentre, iside);
	dj = MAX(jcentre, jside) - MIN(jcentre, jside);
	if (!di && !dj) {
		fprintf(stderr, "doplot: cannot plot zero size selection\n");
		return;
	}
	N = 2 * MAX(ABS(di), ABS(dj));
	io = icentre - N / 2;
	jo = jcentre - N / 2;
	
	switch (plottype) {
		case CONTOURPLOT:
			strcpy(plotcommand, "contour");
			break;
		case PROFILEPLOT:
			sprintf(plotcommand, "profile -c %d %d", N / 2, N / 2);
			break;
		case CUSTOMPLOT:
			strcpy(plotcommand, customcommand);
			break;
	}

	fprintf(stderr, "plotting %d x %d selection io=%d jo=%d with command %s\n",
		N, N, io, jo, plotcommand);
		
	ff = (short *) calloc(N, sizeof(short));
	if (!(fitspipe = popen(plotcommand, "w")))
		error_exit("doplot: failed to open pipe\n");
	fits = new2Dfitsheader(N, N, SHORT_PIXTYPE);
	fits->opstream = fitspipe;
	fits->intpixtype = SHORT_PIXTYPE;
	com = newtextcomment("HISTORY" , "created by xim", NULL);
	appendcomment(com, fits);
	sprintf(astring, "command: %s", plotcommand);
	com = newtextcomment("HISTORY" , astring, NULL);
	appendcomment(com, fits);
	sprintf(astring, "%d x %d selection: io=%d jo=%d", N, N, io, jo);
	com = newtextcomment("HISTORY" , astring, NULL);
	appendcomment(com, fits);
	
	writefitsheader(fits);
	for (i = io; i < io + N; i++) {
		for (j = jo; j < jo + N; j++) {
			if (i < 0 || i >= N2[0] || j < 0 || j >= N1[0])
				ff[j - jo] = MAGIC;
			else
				ff[j - jo] = f[0][0][i][j];
		}
		writefitsline(ff, fits);
	}
	writefitstail(fits);
	pclose(fitspipe);
	free(ff);
}









