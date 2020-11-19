/** Header file generated with fdesign on Sun Mar 24 19:21:57 2002.**/

#ifndef FD_control_panel_h_
#define FD_control_panel_h_

/** Callbacks, globals and object handlers **/
extern void rewind_button_callback(FL_OBJECT *, long);
extern void stop_button_callback(FL_OBJECT *, long);
extern void rec_button_callback(FL_OBJECT *, long);
extern void forwardstep_button_callback(FL_OBJECT *, long);
extern void quit_button_callback(FL_OBJECT *, long);
extern void command_input_callback(FL_OBJECT *, long);
extern void zscale_counter_callback(FL_OBJECT *, long);
extern void go_button_callback(FL_OBJECT *, long);


/**** Forms and Objects ****/
typedef struct {
	FL_FORM *control_panel;
	void *vdata;
	char *cdata;
	long  ldata;
	FL_OBJECT *rewind_button;
	FL_OBJECT *stop_button;
	FL_OBJECT *rec_button;
	FL_OBJECT *forward_step_button;
	FL_OBJECT *frame_number_text;
	FL_OBJECT *quit_button;
	FL_OBJECT *command_input;
	FL_OBJECT *zscale_counter;
	FL_OBJECT *go_button;
} FD_control_panel;

extern FD_control_panel * create_form_control_panel(void);

#endif /* FD_control_panel_h_ */
