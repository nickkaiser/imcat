/*
 * fits3Dviewer.c: geomview module
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>           /* for struct timeval below */
#include <sys/types.h>
#include <unistd.h>

#include "forms.h"              /* for FORMS library */
#define	TRUE 1
#include "xforms-compat.h" 	/* defines foreground - to do nothing! */
#include "imlib/fits.h"
#include "utils/arrays.h"
#include "control_panel.h"
#include "callbacks.h"
#include "mesh.h"

#define usage "\nNAME\n\
        fits3Dviewer - use geomview to display 3D FITS files\n\
\n\
SYNOPSIS\n\
        fits3Dviewer\n\
\n\
DESCRIPTION\n\
        fits3Dviewer is a geomview module.  To use it,\n\
	move to a working directory, and run the script\n\
	geomviewsetup.  Then start geomview.  'FITS 3D viewer'\n\
	should appear in the 'External Modules' menu.  Clicking on\n\
	this should bring up a window 'fits3Dviewer' into which you can enter\n\
	an arbitrary command that, when executed, will generate\n\
	a 3-D FITS file.  The planes of that image will be rendered\n\
	by geomview as a series of surface plots.  This is a useful way\n\
	to make animations.  One can use the controls in the various\n\
	windows to fine-tune the appearance of the surface and the\n\
	z-scale.  Then, pressing the 'record' followed by the 'play'\n\
	button will result in a series of files /tmp/fits3Dviewer.nnnn.ppm\n\
	being created.  These can then be converted into e.g. animated\n\
	GIF format using the gimp.\n\
\n\
	For example\n\
	 ic -c 1024 1024 'grand 3 *' | smooth -p -1 | smooth -f 0 0 2 | helicalscan 64 2 22.5\n\
\n\
SEE ALSO\n\
        geomviewsetup\n\
AUTHOR\n\
        Nick Kaiser --- kaiser@hawaii.edu\n\n"

/* parameters for the wave model */
float			log2zscale, zscale, **f;
int			opmode, recording, framenumber, nframes, Nx, Ny;
FD_control_panel 	*thecontrolpanel;
fitsheader		*fits;
char			*fitsgencommand;
FILE			*ipf;


/* #define ASCII_DATA */

main(argc, argv)        
     char **argv;
{
  /* initialise */
  opmode 		= IDLE_MODE;	/* start idling */
  log2zscale 		= 0.0;		/* scale factor */
  zscale		= 1.0;
  recording 		= 0;
  framenumber 		= 0;
  f			= NULL;
  ipf 			= NULL;

/* check for usage argument  */
if ((argc == 2) && (!strcmp(argv[1], "-u"))) {
	fprintf(stderr, usage);
	exit(1);
}
  /* Forms panel setup.
   */
#ifdef XFORMS
  FL_INITIALIZE("fits3Dviewer");
#else
  foreground();
#endif
  thecontrolpanel = create_form_control_panel();
  fl_set_counter_bounds(thecontrolpanel->zscale_counter, -1000.0, 1000.0);
  fl_set_counter_value(thecontrolpanel->zscale_counter, log2zscale);
  fl_set_counter_filter(thecontrolpanel->zscale_counter, counter_filter);
/* testing
  fl_set_input(thecontrolpanel->command_input, "helicalscan 128 1 45 < /tmp/kolmogorov.fits"); */
  fl_set_focus_object(thecontrolpanel->control_panel, thecontrolpanel->command_input);
  fl_show_form(thecontrolpanel->control_panel, FL_PLACE_SIZE, TRUE, "fts3Dviewer");

  setenv("LD_LIBRARY_PATH", "/usr/local/lib", 1);

  /* Geomview setup.
   */
  printf("(geometry fits3Dviewer { : foo })\n");
  fflush(stdout);

  /* Loop until killed.
   */
  while (1) {	
    fl_check_forms();
    if (opmode == RUN_MODE) {
    	    ReadNewFrame();
	    DrawMesh();
    }
  }
}

