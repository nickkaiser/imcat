#include <stdio.h>
#include <math.h>
#include "forms.h"
#include "xforms-compat.h"
#include "imlib/fits.h"
#include "utils/arrays.h"
#include "control_panel.h"
#include "callbacks.h"
#include "mesh.h"

extern float 		log2zscale, zscale, **f;
extern int		opmode, recording, framenumber, nframes, Nx, Ny, haveimage;
extern FD_control_panel *thecontrolpanel;
extern fitsheader	*fits;
extern char		*fitsgencommand;
extern FILE		*ipf;


/* This is the interface
   between our program and the forms code.
 */

void zscale_counter_callback(FL_OBJECT *obj, long val)
{
  	log2zscale = fl_get_counter_value(thecontrolpanel->zscale_counter);
  	zscale = pow(2.0, log2zscale);
	if (f) {
	  	DrawMesh();
	}
}

void rewind_button_callback(FL_OBJECT *obj, long val)
{
	if (f) {
		start_command();
		ReadNewFrame();
		DrawMesh();
	}
  	opmode = IDLE_MODE;
}


void forwardstep_button_callback(FL_OBJECT *obj, long val)
{
	if (f) {
		ReadNewFrame();
		DrawMesh();
  		opmode = IDLE_MODE;
	}
}

void stop_button_callback(FL_OBJECT *obj, long val)
{
  opmode = IDLE_MODE;
}

void go_button_callback(FL_OBJECT *obj, long val)
{
	if (f) {
  		opmode = RUN_MODE;
	}
}

void rec_button_callback(FL_OBJECT *obj, long val)
{
  recording = (recording ? 0 : 1);
}

void quit_button_callback(FL_OBJECT *obj, long val)
{
  exit(0);
}


void command_input_callback(FL_OBJECT *obj, long val)
{
  	fitsgencommand = (char *) fl_get_input(obj);
	start_command();
	ReadNewFrame();
	DrawMesh();
}


void	start_command(void)
{
	fprintf(stderr, "# starting new command: '%s' \n", fitsgencommand);
	if (ipf) {
		close(ipf);
	}
  	ipf = popen(fitsgencommand, "r");
  	if (!ipf) {
		fprintf(stderr, "# failed to open input pipe using command '%s'\n", fitsgencommand);
		return;
	}
	fits = readfitsheader(ipf);
  	if (!fits) {
		fprintf(stderr, "# failed to read a valid FITS header from command '%s'\n", fitsgencommand);
		return;
	}
	if (fits->ndim != 3) {
		fprintf(stderr, "# FITS header generated by '%s' is not 3 dimensional\n", fitsgencommand);
		return;
	}
	Nx 	= fits->n[0];
	Ny 	= fits->n[1];
	nframes = fits->n[2];
	if (f) {
		freeFloatArray(f, Nx, Ny);
	}
	allocFloatArray(&f, Nx, Ny);
	framenumber = 0;
}


char *counter_filter(FL_OBJECT *obj, double value, int prec)
{
 	static char buf[64];

       /* sprintf(buf, "%.*f", prec, pow(2.0, value)); */
       sprintf(buf, "%8.3lg", pow(2.0, value));
       return buf;
}  
