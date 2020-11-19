#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/extensions/shape.h>
#ifdef IMLIB2
#include <Imlib2.h>
#else
#include <Imlib.h>
#endif
#include        "utils/error.h"
#include        "utils/arrays.h"
#include        "utils/args.h"
#include        "imlib/fits.h"


#define usage "\n\n\
NAME\n\
        xfv --- FITS stream viewer\n\
\n\
SYNOPSIS\n\
        xvf [-b blocksize]\n\
\n\
DESCRIPTION\n\
        xvf is an X-based FITS Viewer for viewing 3D FITS\n\
        images as a movie. It reads from stdin a series of\n\
	8 bit FITS images.  These can either be monochrome,\n\
	in which case the image is 3D: f[NFRAMES][NY][NX], or\n\
	or in RGB format: f[NFRAMES][3][NY][NX] and\n\
	consisting of the R, G, B planes of a set of images.\n\
\n\
	It does no scaling.\n\
\n\
	You can use 'colorize' to convert monochrome to color.\n\
\n\
AUTHOR\n\
        Nick Kaiser --- kaiser@hawaii.edu\n\n"

static unsigned char	*srcdat, *imdat;
int	getimageplane(void);
static	int	w, h, rgbsourceimage;

int main(int argc, char **argv)
{
     	Display *disp;
     	ImlibData *id;
     	XSetWindowAttributes attr;
     	Window win;
     	ImlibImage *im;
     	Pixmap p,m;
 	fitsheader	*fits;
	int		first, nframes, f, blocksize;
     	char		*flag;

	/* defaults */
	blocksize = 1;

	argsinit(argc, argv, usage);
	while (flag = getflag()) {
		switch (flag[0]) {
			case 'b':
				blocksize = getargi();
				break;
			default:
				error_exit(usage);
		}
	}

	fits = readfitsheader(stdin);
	if (fits->extpixtype != UCHAR_PIXTYPE) {
		error_exit("xfv: source image must be 8 bit (unsigned char) format\n");
	}
	switch (fits->ndim) {
		case 4:
			rgbsourceimage = 1;
			nframes = fits->n[3];
			if (fits->n[2] != 3) {
				error_exit("xfv: rgb source image must have M x 3 x N2 x N1 dimensionality\n");
			}
			break;
		case 3:
			rgbsourceimage = 0;
			nframes = fits->n[2];
			break;
		default:
			error_exit("xfv: source image must be 3 or 4 dimensional\n");
	}	
	w = fits->n[0];
	h = fits->n[1];

	if (rgbsourceimage) {
		srcdat = (unsigned char *) calloc(3 * h * w, sizeof(unsigned char));
	} else {
		srcdat = (unsigned char *) calloc(h * w, sizeof(unsigned char));
	}
	imdat = (unsigned char *) calloc(3 * h * w, sizeof(unsigned char));

 	/* Connect to the default Xserver */
    	disp=XOpenDisplay(NULL);
 	/* Immediately afterwards Intitialise Imlib */
    	id=Imlib_init(disp);
 	/* Create a Window to display in */
   	win=XCreateWindow(disp, DefaultRootWindow(disp), 0, 0, blocksize * w, blocksize * h, 0, id->x.depth,
        	InputOutput,id->x.visual, 0, &attr);
    	XSelectInput(disp,win,StructureNotifyMask);

 	first = 1;
	p = (Pixmap) NULL;
	for (f = 0; f < nframes; f++) {
		getimageplane();
		if (first) {
			im = Imlib_create_image_from_data(id, imdat, NULL, w, h);
			/* im = Imlib_create_image_from_data(id, imdat, NULL, h, w); */
		} else {	
			im->rgb_data = imdat;
			Imlib_changed_image(id, im);
		}
		Imlib_render(id, im, blocksize * w, blocksize * h);
		/* Imlib_render(id, im, blocksize * h, blocksize * w); */
		if (p) {
			Imlib_free_pixmap(id,p);
		}
		p=Imlib_move_image(id,im);
		XSetWindowBackgroundPixmap(disp,win,p);
		if (first) {
			XMapWindow(disp,win);
			first = 0;
		} else {
			XClearWindow(disp,win);
		}
		XSync(disp,False);
	}

 /* Event loop to handle resizes */   
    for(;;)
      {
        XEvent ev;
    
 /* Sit and wait for an event to happen */ 
        XNextEvent(disp,&ev);
        if (ev.type==ConfigureNotify)
          {
            w=ev.xconfigure.width;h=ev.xconfigure.height;
 /* Re-render the Image to the new size */ 
            Imlib_render(id,im,w,h);
 /* Free the previous pixmap used for the window - note ONLY the pixmap is */
 /* freed - the mask is marked automatically to be freed along with the */
 /* pixmap. There is no need to free it as well - in fact it is advised you do */
 /* not. You must use the Imlib free function because of Imlib's caching. Do */
 /* not use any other free functions. You can use this function for Pixmaps */
 /* not created by Imlib - and it will just go free them as expected. */
            Imlib_free_pixmap(id,p);
            p=Imlib_move_image(id,im);
 /* The mask will be 0 if the image has no transparency */
            m=Imlib_move_mask(id,im);
 /* Put the Image pixmap in the background of the window */
            XSetWindowBackgroundPixmap(disp,win,p);
 /* If there was a mask to the image, set the Image's mask to it */
            if (m) XShapeCombineMask(disp,win,ShapeBounding,0,0,m,ShapeSet);
 /* Clear the window to update the background change */
            XClearWindow(disp,win);
 /* Synchronise with the Xserver */
            XSync(disp,False);
          }
      }
 }


int	getimageplane(void)
{
	int	c, x, y;

	if (rgbsourceimage) {
		if (!fread(srcdat, sizeof(char), 3 * h * w, stdin)) {
			return(0);
		}
		for (c = 0; c < 3; c++) {
			for (y = 0; y < h; y++) {
				for (x = 0; x < w; x++) {
					imdat[3 * (y * w + x) + c] = srcdat[h * w * c + (h - y - 1) * w + x];
				}
			}
		}
	} else {
		if (!fread(srcdat, sizeof(char), h * w, stdin)) {
			return(0);
		}
		for (y = 0; y < h; y++) {
			for (x = 0; x < w; x++) {
				for (c = 0; c < 3; c++) {
					imdat[3 * (y * w + x) + c] = srcdat[(h - y - 1) * w + x];
				}
			}
		}
	}
	return(1);
}
