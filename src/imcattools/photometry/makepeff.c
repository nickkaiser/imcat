#define	usage "\n\n\n\
NAME\n\
	makepeff --- compute effective polarisability\n\
\n\
SYNOPSIS\n\
	makepeff [options...]\n\
\n\
DESCRIPTION\n\
	'makepeff' first reads from stdin a catalogue which\n\
	must contain at least the following entries:\n\
		F	# windowed flux\n\
		q0	# size\n\
		q[2]	# polarisation\n\
		R[2]	# flux response\n\
		P0[2]	# size response\n\
		P[2][2]	# polarisation response\n\
	as created by 'getshapes2'.\n\
	It bins appropriate combinations of these in a cubical\n\
	array in F, p0, q space and computes P_effective.\n\
\n\
OPTIONS\n\
	Options are\n\
		-u			# print this message\n\
		-F logF1 logF2 nF	# range of log_10 F and number of bins defaults to 2 3.6 8\n\
		-q q1 q2 nq		# range of |q| and number of bins defaults to 0.0 0.5 32\n\
		-Q q01 q02 nq0		# range of q0 and number of bins defaults to 2.5 3.5 32\n\
\n\
OUTPUT\n\
	Output is a 3-plane FITS file with planes containing\n\
		plane 0		n	# number of objects in cell = sum 1\n\
		plane 1		nP	# sum P_eff\n\
		plane 2		nq	# sum sqrt(q.q)\n\
		plane 3		nPbar	# sum P\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@hawaii.edu\n\
\n\n\n"		


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../utils/ipbuff.h"
#include "../../utils/arrays.h"
#include "../../utils/args.h"

float	***alloc3Darray(int N1, int N2, int N3);
float	***boxavg3D(float ***f, int N1, int N2, int N3);
float	***grad3D(float ***f, int dir, double dx, int N1, int N2, int N3);

main(int argc, char *argv[])	
{
  double		logF, logFmin, logFmax, dlogF, q0, q0min, q0max, dq0, q, qmin, qmax, dq;
  double		qq, P, P0q, Rq, qPq;
  int		nF, nq0, nq, iF, iq0, iq;
  double		*buff;
  int		buffwidth, dim[4];
  FILE		*ipf;
  float		***n, ***nP, ***nqPqhat, ***nP0q, ***nRqhat;
  float		***ns, ***nPs, ***gradnqPqhat, ***gradnP0q, ***nPeff, ***nqq, ***nqs;
  float		***gradnRqhat_F, ***gradnRqhat_q0, ***gradnRqhat_q, ***nRqeff;
  fitsheader	*fits;
  char		*flag;

  /* defaults */
  logFmin = 2.0;
  logFmax = 3.6;
  nF = 9;
  q0min = 2.5;
  q0max = 3.5;
  nq0 = 32;
  qmin = 0.0;
  qmax = 0.5;
  nq = 32;
	
 	/* parse args */
	argsinit(argc, argv, usage);
	while(flag = getflag()) {
		switch (flag[0]) {
			case 'F':
				logFmin = (double) getargf();
				logFmax = (double) getargf();
				nF = getargi();
				break;
			case 'q':
				qmin = (double) getargf();
				qmax = (double) getargf();
				nq = getargi();
				break;
			case 'Q':
				q0min = (double) getargf();
				q0max = (double) getargf();
				nq0 = getargi();
				break;
			case 'u':
			default:
				error_exit(usage);
		}
	}
	


  /* compute bin sizes */
  dlogF 	= (logFmax - logFmin) / nF;
  dq0	= (q0max - q0min) / nq0;
  dq	= (qmax - qmin) / nq;

  /* open catalogue */
  ipf = popen("lc -b -o 'logF = %F log10' q0 'qq = %q %q dot' 'P = %P[0][0] %P[1][1] + 2 /' 'P0q = %P0 %q dot 2 /' 'Rq = %R %q dot 2 /' 'qPq = %q %P %q dot dot 2 /'", "r");
  if (!ipf) {
    error_exit("makepeff: failed to open lc input pipe\n");
  }
  buffwidth = 7;
  buff = (double *) calloc(buffwidth, sizeof(double));

  /* allocate arrays */
  n 	= alloc3Darray(nq, nq0, nF);
  nP	= alloc3Darray(nq, nq0, nF);
  nqPqhat	= alloc3Darray(nq, nq0, nF);
  nP0q	= alloc3Darray(nq, nq0, nF);
  nRqhat	= alloc3Darray(nq, nq0, nF);
  nRqeff	= alloc3Darray(nq, nq0, nF);
 nPeff	= alloc3Darray(nq, nq0, nF);
  nqq	= alloc3Darray(nq, nq0, nF);
  nqs	= alloc3Darray(nq, nq0, nF);

  /* accumulate arrays */
  while (fread((void *) buff, sizeof(double), buffwidth, ipf) == buffwidth) {
    logF 	= buff[0];
    q0	= buff[1];
    qq	= buff[2];
    P	= buff[3];
    P0q	= buff[4];
    Rq	= buff[5];
    qPq	= buff[6];
    q = sqrt(qq);
    if (q > FLT_MIN) {
      iF = floor((logF - logFmin) / dlogF);
      iq0 = floor((q0 - q0min) / dq0);
      iq = floor((q - qmin) / dq);
      if ((iF >= 0 && iF < nF) && (iq0 >= 0 && iq0 < nq0) && (iq >= 0 && iq < nq)) {
	n[iF][iq0][iq] 		+= 1;
	nqq[iF][iq0][iq] 	+= q;
	nP[iF][iq0][iq]		+= P;
	nqPqhat[iF][iq0][iq]	+= qPq / q;
	nP0q[iF][iq0][iq]	+= P0q;
      	nRqhat[iF][iq0][iq]	+= Rq / q;
      }
    }
  }

  /* make smoothed arrays */
  ns 	= boxavg3D(n, nq, nq0, nF);
  nPs 	= boxavg3D(nP, nq, nq0, nF);
  nqs 	= boxavg3D(nqq, nq, nq0, nF);

  /* make the gradient fields */
  gradnqPqhat = grad3D(nqPqhat, 1, dq, nq, nq0, nF);
  gradnP0q = grad3D(nP0q, 2, dq0, nq, nq0, nF);
  gradnRqhat_F = grad3D(nRqhat, 3, dlogF, nq, nq0, nF); 
  gradnRqhat_q0 = grad3D(nRqhat, 2, dq0, nq, nq0, nF);
  gradnRqhat_q = grad3D(nRqhat, 1, dq, nq, nq0, nF);

  /* compute the combined Rq term */
  for (iF = 0; iF < nF; iF++) {
    logF = logFmin + (iF + 0.5) * dlogF;
    for (iq0 = 0; iq0 < nq0; iq0++) {
      q0 = q0min + (iq0 + 0.5) * dq0;
      for (iq = 0; iq < nq; iq++) {
	q = qmin + (iq + 0.5) * dq;
	nRqeff[iF][iq0][iq] = q * (2 * nRqhat[iF][iq0][iq]
				   - gradnRqhat_F[iF][iq0][iq] / log(10.0)
				   + q0 * gradnRqhat_q0[iF][iq0][iq]
				   + q * gradnRqhat_q[iF][iq0][iq]);
      }
    }
  }
	
  /* compute nPeff */
  for (iF = 0; iF < nF; iF++) {
    for (iq0 = 0; iq0 < nq0; iq0++) {
      for (iq = 0; iq < nq; iq++) {
	nPeff[iF][iq0][iq] = nPs[iF][iq0][iq] - gradnqPqhat[iF][iq0][iq] 
	  - gradnP0q[iF][iq0][iq] + nRqeff[iF][iq0][iq];
      }
    }
  }
  nPeff = boxavg3D(nPeff, nq, nq0, nF);

  /* write the array */
  dim[0] = nq;
  dim[1] = nq0;
  dim[2] = nF;
  dim[3] = 4;
  fits = newfitsheader(4, dim, FLOAT_PIXTYPE);
  add_comment(argc, argv, fits);
  appendcomment(newnumericcomment("logFmin", (double) logFmin, NULL), fits);
  appendcomment(newnumericcomment("logFmax", (double) logFmax, NULL), fits);
  appendcomment(newnumericcomment("q0min", (double) q0min, NULL), fits);
  appendcomment(newnumericcomment("q0max", (double) q0max, NULL), fits);
  appendcomment(newnumericcomment("qmin", (double) qmin, NULL), fits);
  appendcomment(newnumericcomment("qmax", (double) qmax, NULL), fits);
  writefitsheader(fits);
  writefitscube((void ***) ns, fits);
  writefitscube((void ***) nPeff, fits);
  writefitscube((void ***) nqs, fits);
  writefitscube((void ***) nPs, fits);
  writefitstail(fits);

  /* all done */
  exit(0);
}


float	***alloc3Darray(int N1, int N2, int N3)
{
  int	size, i, j;
  float	***f;

  size = N1 * N2 * N3;
  f = (float ***) calloc(N3, sizeof(float**));
  for (i = 0; i < N3; i++) {
    f[i] = (float **) calloc(N2, sizeof(float));
    for (j = 0; j < N2; j++) {
      if (!i && !j) {
	f[0][0] = (float *) calloc(size, sizeof(float));
	if (!(f[0][0])) {
	  error_exit("alloc3Darray: calloc failed\n");
	}	
      } else {
	f[i][j] = f[0][0] + N1 * (N2 * i + j);
      }
    }
  }
  return(f);	
}

float	***boxavg3D(float ***f, int N1, int N2, int N3)
{
  int	i1, i2, i3, d1, d2, d3;
  float	***fs;

  fs = alloc3Darray(N1, N2, N3);
  for (i3 = 1; i3 < N3 - 1; i3++) {
    for (i2 = 1; i2 < N2 - 1; i2++) {
      for (i1 = 1; i1 < N1 - 1; i1++) {
	for (d3 = -1; d3 <= 1; d3++) {
	  for (d2 = -1; d2 <= 1; d2++) {
	    for (d1 = -1; d1 <= 1; d1++) {
	      fs[i3][i2][i1] += f[i3+d3][i2+d2][i1+d1];	
	    }
	  }
	}
	fs[i3][i2][i1] /= 27.0;
      }
    }
  }

  return(fs);
}

float	***grad3D(float ***f, int dir, double dx, int N1, int N2, int N3)
{
  int	i1, i2, i3, d1, d2, d3;
  float	***df;

  df = alloc3Darray(N1, N2, N3);
  for (i3 = 1; i3 < N3 - 1; i3++) {
    for (i2 = 1; i2 < N2 - 1; i2++) {
      for (i1 = 1; i1 < N1 - 1; i1++) {
	for (d3 = -1; d3 <= 1; d3++) {
	  for (d2 = -1; d2 <= 1; d2++) {
	    for (d1 = -1; d1 <= 1; d1++) {
	      switch (dir) {
	      case 1:
		df[i3][i2][i1] += d1 * f[i3+d3][i2+d2][i1+d1];
		break;	
	      case 2:
 		df[i3][i2][i1] += d2 * f[i3+d3][i2+d2][i1+d1];
		break;	
	      case 3:
		df[i3][i2][i1] += d3 * f[i3+d3][i2+d2][i1+d1];
		break;
	      default:
		error_exit("grad3D: illegal direction\n");
	      }	
	    }
	  }
	}
	df[i3][i2][i1] /= (18.0 * dx);
      }
    }
  }
  return(df);
}
