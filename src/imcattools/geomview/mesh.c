/*
 * functions to write the mesh
 */
#include <math.h>
#include <stdio.h>
#include "forms.h"
#include "imlib/fits.h"
#include "utils/arrays.h"
#include "mesh.h"
#include "control_panel.h"

extern float		zscale, **f;
int			opmode, recording, framenumber, nframes, Nx, Ny;
extern FD_control_panel *thecontrolpanel;
extern fitsheader	*fits;
extern FILE		*ipf;
static char		framestring[256];

void	ReadNewFrame(void)
{
	if (framenumber < nframes - 1) {
  		/* update the frame number string */
  		sprintf(framestring, "frame: %5d", framenumber);
  		fl_set_object_label(thecontrolpanel->frame_number_text, framestring);
		/* read the next plane */
		readfitsplane((void *) f, fits);
		framenumber++;		
	}
}

void	DrawMesh(void)
{
  int 	i,j;
  char	bindat[12], snapshotstring[256];

  printf("(read geometry { define foo \n");
#ifdef ASCII_DATA
  printf("MESH\n");
  printf("%4d %4d\n", Nx, Ny);
#else
  printf("MESH BINARY");
  ((int *) bindat)[0] = Nx;
  ((int *) bindat)[1] = Ny;
  swapbytes((void *) bindat, 2);
  fwrite(bindat, sizeof(char), 8, stdout);
#endif
  for (j = 0; j < Ny; j++) {
    for (i = 0; i < Nx; i++) {
	#ifdef ASCII_DATA
		printf("%f %f %f\t", (float) i, (float) j, zscale * f[j][i]);
	#else
  		((float *) bindat)[0] = (float) i;
  		((float *) bindat)[1] = (float) j;
  		((float *) bindat)[2] = zscale * f[j][i];
  		swapbytes((void *) bindat, 3);
  		fwrite(&bindat, sizeof(char), 12, stdout);
	#endif
    }
    #ifdef ASCII_DATA
	printf("\n");
    #endif
  }
  printf("})\n");
  if (recording) {
    sprintf(snapshotstring, "(snapshot targetcam /tmp/fits3Dviewer.%04d.ppm ppm)\n", framenumber - 1);
    fprintf(stderr, "# writing frame to disk: %s", snapshotstring);
    printf(snapshotstring);
  }
  fflush(stdout);
}

void	swapbytes(void *addr, int nint)
{
	static char	tmpchar[4];
	int		i, b;
	char		*caddr;

	caddr = (char *) addr;
	for (i = 0; i < nint; i++) {
		for (b = 0 ; b < 4; b++) {
			tmpchar[b] = caddr[b];
		}
		for (b = 0 ; b < 4; b++) {
			caddr[b] = tmpchar[3 - b];
		}
		caddr += 4;
	}
}

