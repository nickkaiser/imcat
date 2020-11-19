/*
 * front.c
 *
 * creates the X Fits Viewer application xfv.
 * This interface realised with Athena widgets.
 *
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

#include <X11/Xaw/Label.h>	
#include <X11/Xaw/Box.h>	
#include <X11/Xaw/Form.h>	
#include <X11/Xaw/Command.h>	
#include <X11/Xaw/Viewport.h>	
#include <X11/Xaw/Scrollbar.h>	
#include <X11/Xaw/Dialog.h>	
#include <X11/Xaw/Cardinals.h>	

#include	"front.h"
#include	"back.h"
#include	"../../utils/error.h"
#include	"../../utils/arrays.h"
#include	"../../imlib/fits.h"


/* resource structure to hold app-specific resources */
typedef struct _AppResources {
    float fmin, fmax;
	int	blocksize, colormapindex;
} AppResources;

/* resource specification: name, class, representation type, location in structure, default */
static XtResource resources[] = {
    {"fmin", "Fmin", XtRFloat, sizeof(float),
     XtOffsetOf(AppResources, fmin), XtRString, "0.0" },
    {"fmax", "Fmax", XtRFloat, sizeof(float),
     XtOffsetOf(AppResources, fmax), XtRString, "256.0" },
    {"blocksize", "Blocksize", XtRInt, sizeof(int),
     XtOffsetOf(AppResources, blocksize), XtRImmediate, (XtPointer) 1},
    {"colormapindex", "Colormapindex", XtRInt, sizeof(int),
     XtOffsetOf(AppResources, colormapindex), XtRImmediate, (XtPointer) -1},
};

String fallback_resources[] = { 
 	"*Viewport*allowVert: True",
	"*Viewport*allowHoriz: True",
       	"*Scrollbar*length:	100",
	"*picwidget.cursor:	crosshair",
	"*picwidget.foreground:	yellow",
	"*Form.background:	DarkSlateBlue",
	"*Command.background:	CornflowerBlue",
	"*Label.background:	LightBlue",
	"*Scrollbar.background:	LightBlue",
	"*picwidget.background:	grey",
    NULL,
};

static XrmOptionDescRec options[] = {
    {"-min", "fmin",   XrmoptionSepArg,       NULL},
    {"-max", "fmax",   XrmoptionSepArg,        NULL},
    {"-b", "blocksize",   XrmoptionSepArg,        NULL},
    {"-c", "colormapindex",   XrmoptionSepArg,        NULL},
};

/* help strings */
static char *helpmsg[] = {
    "-min fmin              	minimum f-value (0.0)",
    "-max fmax              	maximum f-value (255)",
    "-b   blocksize             paint pixels this size (1)",
    "-c   cmapindex             colormap (-1)",
NULL };


static	Widget	toplevel;
static 	Widget	picwidget, viewport;
static	Widget	topform, stats;
static	Widget	flabel, tlabel;

static	Widget	widgetlist[100];
static	int	numberlist[100];

/* "local" globals */
static	AppResources	app_resources;
static	char		argstring[256];


main(int argc, char **argv)
{
	XtAppContext 	app_con;
	int		N1, N2;


	argsToString(argc, argv, argstring);
		
	/* create the toplevel widget */
	toplevel = XtAppInitialize(&app_con, "Xim_aw", options, XtNumber(options),
			       &argc, argv, fallback_resources,
			       NULL, ZERO);

	/* XtAppInitialize passes back unused args */
	if (argc != 1) {		
		Syntax(app_con, argv[0]);
	}

	/* load the applications resources */
	XtGetApplicationResources( toplevel, (XtPointer) &app_resources,
			       resources, XtNumber(resources), NULL, ZERO );

	setblocksize(app_resources.blocksize);
	setfrange(app_resources.fmin, app_resources.fmax);
	fprintf(stderr, "fmin fmax %f %f\n", app_resources.fmin, app_resources.fmax);
	setcolormapindex(app_resources.colormapindex);
	readimheader(&N1, &N2);
	
				   
	/* create a form widget to hold controls and viewport */
	topform = XtCreateManagedWidget("topform", formWidgetClass, toplevel,
					    NULL, ZERO);

	/* install the view port */		
	viewport = XtVaCreateManagedWidget("viewport", viewportWidgetClass, topform,
					XtNfromVert, 	NULL,
				    NULL);

	/* create labels */
	stats = XtVaCreateManagedWidget("framelabel", formWidgetClass, topform,
					XtNfromVert, 	viewport,
				    NULL);

	flabel = XtVaCreateManagedWidget("framelabel", labelWidgetClass, stats,
					XtNlabel,	"f=     ",
					XtNfromHoriz, 	NULL,
				    NULL);

	tlabel = XtVaCreateManagedWidget("timelabel", labelWidgetClass, stats,
					XtNlabel,	"t=            ",
					XtNfromHoriz, 	flabel,
				    NULL);


	
	/* now we set up the drawing area: picwidget of commandWidgetClass */
	picwidget = XtVaCreateWidget("picwidget", commandWidgetClass, viewport,
					XtNlabel, 	"",
					XtNwidth,	app_resources.blocksize * N1,
					XtNheight,	app_resources.blocksize * N2,
					XtNresize,	True,
					XtNinternalHeight,	0,
					XtNinternalWidth,	0,
				    NULL);
	XtManageChild(picwidget);
 
	/*
	 * Add a callback routine to the picwidget widget that will print
	 * out the coordinates and pixel value.
	XtAddCallback(picwidget, XtNcallback, picclick, NULL);
	 */

	/* add the background tasks	
	XtAppAddWorkProc(app_con, myworkproc3, NULL);
	XtAppAddWorkProc(app_con, myworkproc2, NULL); */
	XtAppAddWorkProc(app_con, myworkproc1, NULL);
	makewidgetdatabase();
	
	XtRealizeWidget(toplevel);
	alloc_shades();			/* sets many globals like mydisplay etc. */
	makeXImage();			/* make the XImage */

	/* now create the pixmaps */
	create_pixmap();

	/* fill and show the first image */
	fill_pixmap();
	showpixmap();

	/* start the event loop */
	XtAppMainLoop(app_con);
}




void	Syntax(XtAppContext app_con, char *call)
{
	char **msg;

    XtDestroyApplicationContext(app_con);
    fprintf (stderr, "\nNAME\n\txfv --- FITS stream viewer\n\nSYNOPSIS\n\txvf [-options]\n\n");
    fprintf (stderr, "DESCRIPTION\n\txvf is an X-based FITS Viewer for viewing 3D FITS\n");
    fprintf (stderr, "\timages as a movie.\n\n");
    fprintf (stderr, "\tOptions include:\n");
    for (msg = helpmsg; *msg; msg++) fprintf (stderr, "\t\t%s\n", *msg);
    fprintf (stderr, "\nAUTHOR\n\tNick Kaiser --- kaiser@hawaii.edu\n");
   exit(1);
}




void	xfgetargstring(char **argstrptr)
{
	*argstrptr = argstring;
}




float	xfgetvalue(int name)
{
	switch (name) {
		case FMAX:
			return(app_resources.fmax);
			break;
		case FMIN:
			return(app_resources.fmin);
			break;
		default:
			error_exit("xfgetvalue: unknown value name\n");
			break;
	}
}




void	makewidgetdatabase(void)
{
	widgetlist[0] = toplevel;
	numberlist[0] = TOPLEVEL;
 	widgetlist[1] = picwidget;
	numberlist[1] = PICWIDGET;
	widgetlist[2] = flabel;
	numberlist[2] = FLABEL;
	widgetlist[3] = tlabel;
	numberlist[3] = TLABEL;
	widgetlist[4] = NULL;
	numberlist[4] = 0;
}



int	xfgetwidgetnumber(Widget thewidget)
{
	int	index = 0;

	while (1) {
		if (!numberlist[index])
			error_exit("xfgetwidgetnumber: widget not in database\n");
		if (widgetlist[index] == thewidget)
			return(numberlist[index]);
		index++;
	}
}




Widget	xfgetwidget(int thenumber)
{
	int	index = 0;

	while (1) {
		if (!numberlist[index])
			error_exit("xfgetwidget: widget not in database\n");
		if (numberlist[index] == thenumber)
			return(widgetlist[index]);
		index++;
	}
}




void	setlabelstring(Widget w, char *thestring)
{
	Cardinal		n;
	Arg		arg[5];

	n = 0;
	XtSetArg(arg[n], XtNlabel, thestring);	n++;
	XtSetValues(w, arg, n);
}




void	setlabelnumber(Widget w, float thenumber)
{
	char	thestring[20];

	sprintf(thestring, "%.1f\n", thenumber);
	setlabelstring(w, thestring);
}






