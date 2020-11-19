/*
 * zap.h ---- functions to zap/restore image
 */

#define	ZAP_MODE 	0
#define UNZAP_MODE	1

void	zap(int mode, double rmax, double *x, float **f, float **fzap, short **nzap, int N1, int N2);


