/*
 * ximfront.c
 *
 * creates the X-image-viewer application.
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

#include	"ximfront.h"
#include	"ximback.h"
#include	"ximbuttons.h"
#include	"ximcolor.h"
#include	"../../utils/error.h"

/* resource structure to hold app-specific resources */
typedef struct _AppResources {
    int fmin, fmax;
    int goodwindowsize;
} AppResources;

/* resource specification: name, class, representation type, location in structure, default */
static XtResource resources[] = {
    {"fmin", "Fmin", XtRInt, sizeof(int),
     XtOffsetOf(AppResources, fmin), XtRImmediate, (XtPointer) 0},
    {"fmax", "Fmax", XtRInt, sizeof(int),
     XtOffsetOf(AppResources, fmax), XtRImmediate, (XtPointer) 256},
    {"goodwindowsize", "Goodwindowsize", XtRInt, sizeof(int),
     XtOffsetOf(AppResources, goodwindowsize), XtRImmediate, (XtPointer) 512},
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
    {"-goodsize", "goodwindowsize",   XrmoptionSepArg,        NULL},
    {"-min", "fmin",   XrmoptionSepArg,       NULL},
    {"-max", "fmax",   XrmoptionSepArg,        NULL},
};

/* help strings */
static char *helpmsg[] = {
    "-min fmin              	minimum f-value (0)",
    "-max fmax              	maximum f-value (256)",
    "-goodsize                  good window size (512)",
NULL };


/* some parameters */
#define	MAX_SLIDERS	5

static	Widget	toplevel;
static 	Widget	picwidget, viewport;
static	Widget	topform, title, controls, stats, buttons, sliders;
static	Widget	slider[MAX_SLIDERS], scrollbar[MAX_SLIDERS];
static	Widget	sliderlabel[MAX_SLIDERS], slidervalue[MAX_SLIDERS];
static	Widget	zoom, unzoom, makemask, print, plot;
static	Widget	poslabel, pixvallabel, pixrangelabel;
static	Pixmap	pixmap[MAX_LEVELS];

static	Widget	widgetlist[100];
static	int	numberlist[100];

/* "local" globals */
static	AppResources	app_resources;
static	char		argstring[256];


main(int argc, char **argv)
{
	XtAppContext 	app_con;
	int	color;

	argsToString(argc, argv, argstring);
		
	/* create the toplevel widget */
	toplevel = XtAppInitialize(&app_con, "Xim_aw", options, XtNumber(options),
			       &argc, argv, fallback_resources,
			       NULL, ZERO);

	/* XtAppInitialize passes back unused args */
	if (argc != 1)		
		Syntax(app_con, argv[0]);
	
	/* load the applications resources */
	XtGetApplicationResources( toplevel, (XtPointer) &app_resources,
			       resources, XtNumber(resources), NULL, ZERO );
				   
	/* create a form widget to hold controls and viewport */
	topform = XtCreateManagedWidget("topform", formWidgetClass, toplevel,
					    NULL, ZERO);

	/* create contols form to contain title, stats, buttons, sliders */
	controls = XtVaCreateManagedWidget("controls", formWidgetClass, topform,
					XtNfromHoriz, 	NULL,
				    NULL);

	/* top entry in controls is the application title */
	title = XtVaCreateManagedWidget("title", commandWidgetClass, controls,
					XtNlabel,	"xim        ",
					XtNfromVert, 	NULL,
				    NULL);
	
	/* if title is pressed, we pop up the credits box */
	XtAddCallback(title, XtNcallback, Title, NULL);
	

	/* put stats window below the title */					
	stats = XtVaCreateManagedWidget("stats", formWidgetClass, controls,
					XtNfromVert, 	title,
				    NULL);

	poslabel = XtVaCreateManagedWidget("poslabel", labelWidgetClass, stats,
					XtNlabel,	"pos:          ",
					XtNfromVert, 	NULL,
				    NULL);

	pixvallabel = XtVaCreateManagedWidget("pixvallabel", labelWidgetClass, stats,
					XtNlabel,	"pixval:          ",
					XtNfromVert, 	poslabel,
				    NULL);

	/* put button window below the stats */					
	buttons = XtVaCreateManagedWidget("buttons", formWidgetClass, controls,
					XtNfromVert, 	stats,
				    NULL);

	zoom = XtVaCreateManagedWidget("zoom", commandWidgetClass, buttons,
					XtNlabel,	"zoom",
					XtNfromVert, 	NULL,
					XtNfromHoriz, 	NULL,
				    NULL);

	XtAddCallback(zoom, XtNcallback, Zoom, NULL);

	unzoom = XtVaCreateManagedWidget("unzoom", commandWidgetClass, buttons,
					XtNlabel,	"unzoom",
					XtNfromVert, 	NULL,
					XtNfromHoriz, 	zoom,
				    NULL);

	XtAddCallback(unzoom, XtNcallback, UnZoom, NULL);

	makemask = XtVaCreateManagedWidget("makemask", commandWidgetClass, buttons,
					XtNlabel,	"makemask...",
					XtNfromVert, 	zoom,
					XtNfromHoriz, 	NULL,
				    NULL);
					
	XtAddCallback(makemask, XtNcallback, makepopup, NULL);

	print = XtVaCreateManagedWidget("print", commandWidgetClass, buttons,
					XtNlabel,	"postscript...",
					XtNfromVert, 	makemask,
					XtNfromHoriz, 	NULL,
				    NULL);
					
	XtAddCallback(print, XtNcallback, makepopup, NULL);

	plot = XtVaCreateManagedWidget("plot", commandWidgetClass, buttons,
					XtNlabel,	"plot...",
					XtNfromVert, 	print,
					XtNfromHoriz, 	NULL,
				    NULL);
					
	XtAddCallback(plot, XtNcallback, makepopup, NULL);

	/* put sliders window below the stats */					
	sliders = XtVaCreateManagedWidget("sliders", formWidgetClass, controls,
					XtNfromVert, 	buttons,
				    NULL);
	makeslider(0, NULL, NULL);
	makeslider(1, NULL, slider[0]);
	makeslider(2, NULL, slider[1]);

	/* finally we install the view port to the right of the controls */		
	viewport = XtVaCreateManagedWidget("viewport", viewportWidgetClass, topform,
					XtNfromHoriz, 	controls,
				    NULL);

	
	/* now we set up the drawing area: picwidget of commandWidgetClass */
	picwidget = XtVaCreateWidget("picwidget", commandWidgetClass, viewport,
					XtNlabel, 	"",
					XtNwidth,	app_resources.goodwindowsize,
					XtNheight,	app_resources.goodwindowsize,
					XtNresize,	True,
					XtNinternalHeight,	0,
					XtNinternalWidth,	0,
				    NULL);
	XtManageChild(picwidget);
 
	/*
	 * Add a callback routine to the picwidget widget that will print
	 * out the coordinates and pixel value.
	 */
	XtAddCallback(picwidget, XtNcallback, picclick, NULL);

	/* add the background tasks */	
	XtAppAddWorkProc(app_con, myworkproc3, NULL);
	XtAppAddWorkProc(app_con, myworkproc2, NULL);
	XtAppAddWorkProc(app_con, myworkproc1, NULL);
	makewidgetdatabase();
	
	XtRealizeWidget(toplevel);
	XtAppMainLoop(app_con);
}




void	Syntax(XtAppContext app_con, char *call)
{
	char **msg;

    XtDestroyApplicationContext(app_con);
    fprintf (stderr, "\nNAME\n\t%s --- colour image display\n\nSYNOPSIS\n\t%s [-options]\n\n", call, call);
    fprintf (stderr, "DESCRIPTION\n\t%s is an X-based image displayer\n", call);
    fprintf (stderr, "\tIt can display monchrome or\n\tcolour fits images created with makecolourimage\n");
    fprintf (stderr, "\tbut these must be in 16 bit format.\n");
    fprintf (stderr, "\tOptions include:\n");
    for (msg = helpmsg; *msg; msg++) fprintf (stderr, "\t\t%s\n", *msg);
    fprintf (stderr, "\nAUTHOR\n\tNick Kaiser --- kaiser@cita.utoronto.ca\n");
   exit(1);
}




void	makeslider(int index, Widget fromVert, Widget fromHoriz)
{
	/* each slider has a label and a scrollbar*/
	slider[index] = XtVaCreateManagedWidget("slider", formWidgetClass, sliders,
					XtNfromVert, 	fromVert,
					XtNfromHoriz, 	fromHoriz,
				    NULL);
					
	sliderlabel[index] = XtVaCreateManagedWidget("sliderlabel", labelWidgetClass, slider[index],
					XtNlabel,	"   ",
					XtNfromVert, 	NULL,
				    NULL);
	
	scrollbar[index] = XtVaCreateManagedWidget("scrollbar", scrollbarWidgetClass, slider[index],
					XtNfromVert, 	sliderlabel[index],
				    NULL);

 	slidervalue[index] = XtVaCreateManagedWidget("slidervalue", labelWidgetClass, slider[index],
					XtNlabel,	"   ",
					XtNfromVert, 	scrollbar[index],
				    NULL);

	XtAddCallback(scrollbar[index], XtNjumpProc, sliderChanged, NULL);
}




void	xfgetargstring(char **argstrptr)
{
	*argstrptr = argstring;
}




int	xfgetvalue(int name)
{
	switch (name) {
		case FMAX:
			return(app_resources.fmax);
			break;
		case FMIN:
			return(app_resources.fmin);
			break;
		case GOODSIZE:
			return(app_resources.goodwindowsize);
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
	widgetlist[2] = title;
	numberlist[2] = TITLE;
	widgetlist[3] = print;
	numberlist[3] = PRINT;
	widgetlist[4] = makemask;
	numberlist[4] = MAKEMASK;
	widgetlist[5] = poslabel;
	numberlist[5] = POSLABEL;
	widgetlist[6] = pixvallabel;
	numberlist[6] = PIXVALLABEL;
	widgetlist[7] = slider[2];
	numberlist[7] = SLIDER2;
	widgetlist[8] = plot;
	numberlist[8] = PLOT;
	widgetlist[9] = NULL;
	numberlist[9] = 0;
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




void	sliderChanged(Widget w, XtPointer closure, XtPointer call_data)
{
	int	index;
	float	thumbpos;

	thumbpos = *((float *) call_data);
	for (index = 0; index < MAX_SLIDERS; index ++)
		if (scrollbar[index] == w)
			break;
	setsliderval(index, 1 - thumbpos);
	setlabelnumber(slidervalue[index], newslidervalue(index));
   	set_shades();
}





void	setslider(int index, float thumbpos, char *label, char *value)
{
	Arg 		arg[5];
	Cardinal	n;

	n = 0;
	XtSetArg(arg[n], XtNtopOfThumb, (&thumbpos)); n++;
	XtSetValues(scrollbar[index], arg, n);
	setlabelstring(sliderlabel[index], label);
	setlabelstring(slidervalue[index], value);
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






