/*
 * overscancorr.c --- written by Gordon Squires
 */

#define usage "\n\n\n\
NAME\n\
	overscancorr - dark subtraction using overscan region\n\
\n\
SYNOPSIS\n\
	overscancorr [options....] \n\
		-x xmin xmax	# overscan region in x (2052,2086)\n\
		-y ymin ymax	# overscan region in y (0,2047)\n\
\n\
DESCRIPTION\n\
	\"overscancorr\" corrects for the bias dc level. Reads the\n\
	overscan region and fits a linear ramp model along the\n\
	long (y) axis. The dc and gradient is subtracted from the\n\
	image.\n\
\n\
AUTHOR\n\
	Gordon Squires\n\
\n\n\n"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../imlib/fits.h"
#include "../utils/error.h"

#define MAGIC FLOAT_MAGIC

int main(int argc, char *argv[])	
{
  int	   arg = 1, i, j, N1, N2;
  int	   IMIN, IMAX, JMIN, JMAX, pixcount;
  fitsheader	*fits;
  float	   **f, **fcorr;
  float    *rowaverage, a11, a21, a12, a22, detA, b1, b2, b, m, result;

  /* defaults */
  IMIN = 0;
  IMAX = 2047;
  JMIN = 2052;
  JMAX = 2086;
  
  while (arg < argc) {
    if (argv[arg][0] == '-') {
      switch (argv[arg++][1]) {
      case 'y':
	if (1 != sscanf(argv[arg++], "%d", &IMIN))
	  error_exit(usage);
	if (1 != sscanf(argv[arg++], "%d", &IMAX))
	  error_exit(usage);
	break;
      case 'x':
	if (1 != sscanf(argv[arg++], "%d", &JMIN))
	  error_exit(usage);
	if (1 != sscanf(argv[arg++], "%d", &JMAX))
	  error_exit(usage);
	break;
      default:
	error_exit(usage);
	break;
      }
    }
  }
	
  read2Dfloatimage(&f, &N1, &N2, &fits, stdin);
  add_comment(argc, argv, fits); 

  rowaverage = (float *) calloc( (IMAX-IMIN), sizeof(float));

  for (i = IMIN; i <= IMAX; i++) {
    pixcount = 0;
    rowaverage[i-IMIN] = 0.0;
    for (j = JMIN; j <= JMAX; j++) {
      if (f[i][j] != MAGIC)
	{
	  pixcount += 1;
	  rowaverage[i-IMIN] += f[i][j];
	}  
    }
    rowaverage[i-IMIN] /= (float) pixcount;
  }

  a11 = a21 = a12 = a22 = b1 = b2 = 0.0;

  /* calculate matrix elements for max likelihood fit to line */
  for (i = IMIN; i <= IMAX; i++) {  
    a11 += (float) i*i;
    a12 += (float) i;
    b1  += (float) i * rowaverage[i-IMIN];
    b2  += rowaverage[i-IMIN];
  }
  a21 = a12;
  a22 = IMAX - IMIN;

  detA = a11*a22 - a21*a12;
  m    = ( a22*b1  - a12*b2 ) / detA;
  b    = ( -a21*b1 + a11*b2 ) / detA;

  fprintf(stderr,"Gradient fit: f = %f +  %f * i\n", b, m);
  
  for( i = 0; i < N2; i = i + 1 )
    for( j = 0; j < N1; j = j + 1 )
      {
	if( f[i][j] == MAGIC )
	  result = MAGIC;
	else
	  result = f[i][j] - b - m*i; 
	  
	f[i][j] = result;
      }
   
 
  write2Dfloatimage(f, fits); 
  exit(0);
}



