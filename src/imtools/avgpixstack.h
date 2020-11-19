/*
 * avgpixstack.h
 *
 * declaration of avgpixstack() function
 *
 * Nick Kaiser
 *
 * 	avgpixstack() averages a stack of pixels
 *
 *	it reads a vector of pixel values f[nplanes] and associated weight[nplanes]
 *	and averages those with non-zero weight to produce favg and weightsum, the exposure map value
 *	ftemp should be of lenght >= favg
 *
 *	The type of averaging is determined by the opmode parameters
 *
 */

#define	APS_STRAIGHT_AVERAGE	0
#define APS_WEIGHTED_AVERAGE	1
#define APS_MEDIAN		2
#define APS_RANK		3

int	avgpixstack(int opmode, int rank, int nplanes, float *f, float *weight, float *ftemp, float *favg, float *weightsum);
