#define	usage "\n\n\n\
NAME\n\
		imview -- displays a square image\n\
\n\
SYNOPSIS\n\
		imview	[option...] <inputfile\n\
			-f fmin fmax	# max and min f values\n\
			-a 		# ascii file\n\
\n\
DESCRIPTION\n\
		The range of values may be specified in integer units\n\
		otherwise the range is 0 (=white) to 255 (=black)\n\
\n\n\n"		

#include 	<Xm/MainW.h>
#include 	<Xm/Label.h>
#include 	<Xm/Scale.h>
#include 	<Xm/Form.h>
#include 	<stdio.h>
#include	<stdlib.h>
#include	<limits.h>

#include	"imview_stuff.h"
#include	"../../imlib/fits.h"
#include	"../../utils/error.h"

Widget		main_w, fminscale, fmaxscale;
int		fmin = 0, fmax = 256;		
int		n_cells, ascii = 0;

main(int argc, char **argv)
{
	Widget 		toplevel, scrollw, label, form;
	XtAppContext 	app;
	Pixmap		mypixmap;
	short 		**f;
	float		ff;
	int		i, j, arg = 1;	
	char		string[COM_LENGTH * MAX_COMMENTS], argstring[512];

	void	new_fminmax(Widget scale_w, caddr_t client_data, XmScaleCallbackStruct *cbs);

	while (arg < argc) {
		if (*argv[arg] != '-')
			error_exit(usage);
		switch (*(argv[arg++]+1)) {
			case 'f':
				if (!sscanf(argv[arg++], "%d", &fmin))
					error_exit(usage);
				if (!sscanf(argv[arg++], "%d", &fmax))
					error_exit(usage);
				break;
			case 'a':
				ascii = 1;
				break;
			default:
				error_exit(usage);
				break;
		}
	}

		
	/* create Widget interface */		
	toplevel = XtVaAppInitialize(&app, "Demos",
				NULL, 0, &argc, argv, NULL, NULL);

	form = XtVaCreateManagedWidget("form", xmFormWidgetClass, toplevel,
			XmNfractionBase, 8,
			XmNwidth,	500,
			XmNheight,	500,
			NULL);

	fminscale = XtVaCreateManagedWidget("fmin",
		xmScaleWidgetClass, form,
		XtVaTypedArg, 	XmNtitleString, XmRString, "fmin", 4,
		XmNmaximum, 	fmax,
		XmNminimum, 	fmin,
		XmNvalue,	fmin,
		XmNshowValue, 	True,
		XmNtopAttachment,	XmATTACH_POSITION,
		XmNtopPosition,		0,
		XmNleftAttachment,	XmATTACH_POSITION,
		XmNleftPosition,	0,
		XmNrightAttachment,	XmATTACH_POSITION,
		XmNrightPosition,	1,
		XmNbottomAttachment,	XmATTACH_POSITION,
		XmNbottomPosition,	8,
		NULL);

	XtAddCallback(fminscale, XmNdragCallback, new_fminmax, NULL);
	XtAddCallback(fminscale, XmNvalueChangedCallback, new_fminmax, NULL);

	fmaxscale = XtVaCreateManagedWidget("fmax",
		xmScaleWidgetClass, form,
		XtVaTypedArg, 	XmNtitleString, XmRString, "fmax", 4,
		XmNmaximum, 	fmax,
		XmNminimum, 	fmin,
		XmNvalue,	fmax,
		XmNshowValue, 	True,
		XmNtopAttachment,	XmATTACH_POSITION,
		XmNtopPosition,		0,
		XmNleftAttachment,	XmATTACH_POSITION,
		XmNleftPosition,	1,
		XmNrightAttachment,	XmATTACH_POSITION,
		XmNrightPosition,	2,
		XmNbottomAttachment,	XmATTACH_POSITION,
		XmNbottomPosition,	8,
		NULL);

	XtAddCallback(fmaxscale, XmNdragCallback, new_fminmax, NULL);
	XtAddCallback(fmaxscale, XmNvalueChangedCallback, new_fminmax, NULL);

	main_w = XtVaCreateManagedWidget("scrolled_w",
		xmScrolledWindowWidgetClass,   form,
		XmNscrollBarDisplayPolicy, 	XmAS_NEEDED,
        	XmNscrollingPolicy,        	XmAUTOMATIC,
		XmNtopAttachment,	XmATTACH_POSITION,
		XmNtopPosition,		0,
		XmNleftAttachment,	XmATTACH_POSITION,
		XmNleftPosition,	2,
		XmNrightAttachment,	XmATTACH_POSITION,
		XmNrightPosition,	8,
		XmNbottomAttachment,	XmATTACH_POSITION,
		XmNbottomPosition,	8,
       	NULL);


	alloc_grays(&n_cells);
	set_grays(n_cells, 0, n_cells);

	mypixmap = make_pixmap(fmin, fmax, n_cells); 

	label = XtVaCreateManagedWidget("label", xmLabelWidgetClass, main_w,
		XmNlabelType,   XmPIXMAP,
	        XmNlabelPixmap, mypixmap,
	        NULL);

	    /* set the label as the "work area" of the main window */
 	XtVaSetValues(main_w,
	        XmNworkWindow, label,
	        NULL);
	XtRealizeWidget(toplevel);
	XtAppMainLoop(app);
}



void	new_fminmax(Widget scale_w, caddr_t client_data, XmScaleCallbackStruct *cbs)
{
	int	new_fmin, new_fmax, cell_min, cell_max;

	XmScaleGetValue(fminscale, &new_fmin);	
	XmScaleGetValue(fmaxscale, &new_fmax);

	cell_min = (n_cells * (float) (new_fmin - fmin)) / (fmax - fmin);
	cell_max = (n_cells * (float) (new_fmax - fmin)) / (fmax - fmin);
	cell_min = (cell_min < 0 ? 0 : cell_min);
	cell_max = (cell_max > n_cells - 1 ? n_cells - 1 : cell_max);
 	set_grays(n_cells, cell_min, cell_max);	
}



