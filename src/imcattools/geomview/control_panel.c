/* Form definition file generated with fdesign. */

#include "forms.h"
#include <stdlib.h>
#include "control_panel.h"

FD_control_panel *create_form_control_panel(void)
{
  FL_OBJECT *obj;
  FD_control_panel *fdui = (FD_control_panel *) fl_calloc(1, sizeof(*fdui));

  fdui->control_panel = fl_bgn_form(FL_NO_BOX, 571, 201);
  obj = fl_add_box(FL_UP_BOX,0,0,571,201,"");
  fdui->rewind_button = obj = fl_add_button(FL_NORMAL_BUTTON,100,140,30,30,"@4->|");
    fl_set_object_color(obj,FL_DARKER_COL1,FL_DODGERBLUE);
    fl_set_object_lcolor(obj,FL_DODGERBLUE);
    fl_set_object_callback(obj,rewind_button_callback,0);
  fdui->stop_button = obj = fl_add_button(FL_NORMAL_BUTTON,140,140,30,30,"@square");
    fl_set_object_lcolor(obj,FL_RED);
    fl_set_object_callback(obj,stop_button_callback,0);
  fdui->rec_button = obj = fl_add_lightbutton(FL_PUSH_BUTTON,40,140,50,30,"REC");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_callback(obj,rec_button_callback,0);
  fdui->forward_step_button = obj = fl_add_button(FL_NORMAL_BUTTON,180,140,30,30,"@>");
    fl_set_object_color(obj,FL_DARKER_COL1,FL_COL1);
    fl_set_object_lcolor(obj,FL_DODGERBLUE);
    fl_set_object_callback(obj,forwardstep_button_callback,0);
  fdui->frame_number_text = obj = fl_add_text(FL_NORMAL_TEXT,40,90,100,30,"frame:");
    fl_set_object_boxtype(obj,FL_FRAME_BOX);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lalign(obj,FL_ALIGN_LEFT|FL_ALIGN_INSIDE);
  fdui->quit_button = obj = fl_add_button(FL_NORMAL_BUTTON,460,90,70,30,"QUIT");
    fl_set_object_callback(obj,quit_button_callback,0);
  fdui->command_input = obj = fl_add_input(FL_NORMAL_INPUT,40,30,490,40,"command:");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lalign(obj,FL_ALIGN_TOP_LEFT);
    fl_set_object_callback(obj,command_input_callback,0);
  fdui->zscale_counter = obj = fl_add_counter(FL_NORMAL_COUNTER,340,140,190,30,"zscale");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lalign(obj,FL_ALIGN_TOP);
    fl_set_object_callback(obj,zscale_counter_callback,0);
    fl_set_counter_value(obj, 100.0);
    fl_set_counter_return(obj, FL_RETURN_CHANGED);
  fdui->go_button = obj = fl_add_button(FL_NORMAL_BUTTON,220,140,30,30,"@->");
    fl_set_object_color(obj,FL_DARKER_COL1,FL_COL1);
    fl_set_object_lcolor(obj,FL_DODGERBLUE);
    fl_set_object_callback(obj,go_button_callback,0);
  fl_end_form();

  fdui->control_panel->fdui = fdui;

  return fdui;
}
/*---------------------------------------*/

