/*
 * convertarray.c -- written by makeconvertarray.pl
 */

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "fits.h"
#include "convertarray.h"

#define clip(min,max,f) ((f) > (min) ? ((f) < (max) ? (f) : (max)) : (min))

int	convertarray(char *fsrc, char *fdst, int srcpixtype, int dstpixtype, int nel, int bscaling, double bscale, double bzero)
{
	int	i, srcsize, dstsize;

	switch(srcpixtype) {
		case UCHAR_PIXTYPE:
			srcsize = sizeof(unsigned char);
			break;
		case SHORT_PIXTYPE:
			srcsize = sizeof(short);
			break;
		case INT_PIXTYPE:
			srcsize = sizeof(int);
			break;
		case FLOAT_PIXTYPE:
			srcsize = sizeof(float);
			break;
		case DBL_PIXTYPE:
			srcsize = sizeof(double);
			break;
		default:
			error_exit("convertarray: bad pixtype\n");
	}

	switch(dstpixtype) {
		case UCHAR_PIXTYPE:
			dstsize = sizeof(unsigned char);
			switch(srcpixtype) {
				case UCHAR_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((unsigned char *) (fsrc + i * srcsize)) == UCHAR_MAGIC) {
							*((unsigned char *) (fdst + i * dstsize)) = UCHAR_MAGIC;
						} else {
							if (bscaling) {
								*((unsigned char *) (fdst + i * dstsize)) = (unsigned char) clip(0, UCHAR_MAX - 1, floor(0.5 + *((unsigned char *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((unsigned char *) (fdst + i * dstsize)) = (unsigned char) clip(0, UCHAR_MAX - 1, floor(0.5 + *((unsigned char *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case SHORT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((short *) (fsrc + i * srcsize)) == SHORT_MAGIC) {
							*((unsigned char *) (fdst + i * dstsize)) = UCHAR_MAGIC;
						} else {
							if (bscaling) {
								*((unsigned char *) (fdst + i * dstsize)) = (unsigned char) clip(0, UCHAR_MAX - 1, floor(0.5 + *((short *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((unsigned char *) (fdst + i * dstsize)) = (unsigned char) clip(0, UCHAR_MAX - 1, floor(0.5 + *((short *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case INT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((int *) (fsrc + i * srcsize)) == INT_MAGIC) {
							*((unsigned char *) (fdst + i * dstsize)) = UCHAR_MAGIC;
						} else {
							if (bscaling) {
								*((unsigned char *) (fdst + i * dstsize)) = (unsigned char) clip(0, UCHAR_MAX - 1, floor(0.5 + *((int *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((unsigned char *) (fdst + i * dstsize)) = (unsigned char) clip(0, UCHAR_MAX - 1, floor(0.5 + *((int *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case FLOAT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((float *) (fsrc + i * srcsize)) == FLOAT_MAGIC) {
							*((unsigned char *) (fdst + i * dstsize)) = UCHAR_MAGIC;
						} else {
							if (bscaling) {
								*((unsigned char *) (fdst + i * dstsize)) = (unsigned char) clip(0, UCHAR_MAX - 1, floor(0.5 + *((float *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((unsigned char *) (fdst + i * dstsize)) = (unsigned char) clip(0, UCHAR_MAX - 1, floor(0.5 + *((float *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case DBL_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((double *) (fsrc + i * srcsize)) == DBL_MAGIC) {
							*((unsigned char *) (fdst + i * dstsize)) = UCHAR_MAGIC;
						} else {
							if (bscaling) {
								*((unsigned char *) (fdst + i * dstsize)) = (unsigned char) clip(0, UCHAR_MAX - 1, floor(0.5 + *((double *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((unsigned char *) (fdst + i * dstsize)) = (unsigned char) clip(0, UCHAR_MAX - 1, floor(0.5 + *((double *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				default:
					error_exit("convertarray: bad pixtype\n");
			}
			break;
		case SHORT_PIXTYPE:
			dstsize = sizeof(short);
			switch(srcpixtype) {
				case UCHAR_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((unsigned char *) (fsrc + i * srcsize)) == UCHAR_MAGIC) {
							*((short *) (fdst + i * dstsize)) = SHORT_MAGIC;
						} else {
							if (bscaling) {
								*((short *) (fdst + i * dstsize)) = (short) clip(SHRT_MIN + 1, SHRT_MAX, floor(0.5 + *((unsigned char *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((short *) (fdst + i * dstsize)) = (short) clip(SHRT_MIN + 1, SHRT_MAX, floor(0.5 + *((unsigned char *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case SHORT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((short *) (fsrc + i * srcsize)) == SHORT_MAGIC) {
							*((short *) (fdst + i * dstsize)) = SHORT_MAGIC;
						} else {
							if (bscaling) {
								*((short *) (fdst + i * dstsize)) = (short) clip(SHRT_MIN + 1, SHRT_MAX, floor(0.5 + *((short *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((short *) (fdst + i * dstsize)) = (short) clip(SHRT_MIN + 1, SHRT_MAX, floor(0.5 + *((short *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case INT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((int *) (fsrc + i * srcsize)) == INT_MAGIC) {
							*((short *) (fdst + i * dstsize)) = SHORT_MAGIC;
						} else {
							if (bscaling) {
								*((short *) (fdst + i * dstsize)) = (short) clip(SHRT_MIN + 1, SHRT_MAX, floor(0.5 + *((int *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((short *) (fdst + i * dstsize)) = (short) clip(SHRT_MIN + 1, SHRT_MAX, floor(0.5 + *((int *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case FLOAT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((float *) (fsrc + i * srcsize)) == FLOAT_MAGIC) {
							*((short *) (fdst + i * dstsize)) = SHORT_MAGIC;
						} else {
							if (bscaling) {
								*((short *) (fdst + i * dstsize)) = (short) clip(SHRT_MIN + 1, SHRT_MAX, floor(0.5 + *((float *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((short *) (fdst + i * dstsize)) = (short) clip(SHRT_MIN + 1, SHRT_MAX, floor(0.5 + *((float *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case DBL_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((double *) (fsrc + i * srcsize)) == DBL_MAGIC) {
							*((short *) (fdst + i * dstsize)) = SHORT_MAGIC;
						} else {
							if (bscaling) {
								*((short *) (fdst + i * dstsize)) = (short) clip(SHRT_MIN + 1, SHRT_MAX, floor(0.5 + *((double *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((short *) (fdst + i * dstsize)) = (short) clip(SHRT_MIN + 1, SHRT_MAX, floor(0.5 + *((double *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				default:
					error_exit("convertarray: bad pixtype\n");
			}
			break;
		case INT_PIXTYPE:
			dstsize = sizeof(int);
			switch(srcpixtype) {
				case UCHAR_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((unsigned char *) (fsrc + i * srcsize)) == UCHAR_MAGIC) {
							*((int *) (fdst + i * dstsize)) = INT_MAGIC;
						} else {
							if (bscaling) {
								*((int *) (fdst + i * dstsize)) = (int) clip(INT_MIN + 1, INT_MAX, floor(0.5 + *((unsigned char *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((int *) (fdst + i * dstsize)) = (int) clip(INT_MIN + 1, INT_MAX, floor(0.5 + *((unsigned char *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case SHORT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((short *) (fsrc + i * srcsize)) == SHORT_MAGIC) {
							*((int *) (fdst + i * dstsize)) = INT_MAGIC;
						} else {
							if (bscaling) {
								*((int *) (fdst + i * dstsize)) = (int) clip(INT_MIN + 1, INT_MAX, floor(0.5 + *((short *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((int *) (fdst + i * dstsize)) = (int) clip(INT_MIN + 1, INT_MAX, floor(0.5 + *((short *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case INT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((int *) (fsrc + i * srcsize)) == INT_MAGIC) {
							*((int *) (fdst + i * dstsize)) = INT_MAGIC;
						} else {
							if (bscaling) {
								*((int *) (fdst + i * dstsize)) = (int) clip(INT_MIN + 1, INT_MAX, floor(0.5 + *((int *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((int *) (fdst + i * dstsize)) = (int) clip(INT_MIN + 1, INT_MAX, floor(0.5 + *((int *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case FLOAT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((float *) (fsrc + i * srcsize)) == FLOAT_MAGIC) {
							*((int *) (fdst + i * dstsize)) = INT_MAGIC;
						} else {
							if (bscaling) {
								*((int *) (fdst + i * dstsize)) = (int) clip(INT_MIN + 1, INT_MAX, floor(0.5 + *((float *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((int *) (fdst + i * dstsize)) = (int) clip(INT_MIN + 1, INT_MAX, floor(0.5 + *((float *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case DBL_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((double *) (fsrc + i * srcsize)) == DBL_MAGIC) {
							*((int *) (fdst + i * dstsize)) = INT_MAGIC;
						} else {
							if (bscaling) {
								*((int *) (fdst + i * dstsize)) = (int) clip(INT_MIN + 1, INT_MAX, floor(0.5 + *((double *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((int *) (fdst + i * dstsize)) = (int) clip(INT_MIN + 1, INT_MAX, floor(0.5 + *((double *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				default:
					error_exit("convertarray: bad pixtype\n");
			}
			break;
		case FLOAT_PIXTYPE:
			dstsize = sizeof(float);
			switch(srcpixtype) {
				case UCHAR_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((unsigned char *) (fsrc + i * srcsize)) == UCHAR_MAGIC) {
							*((float *) (fdst + i * dstsize)) = FLOAT_MAGIC;
						} else {
							if (bscaling) {
								*((float *) (fdst + i * dstsize)) = (float) (( *((unsigned char *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((float *) (fdst + i * dstsize)) = (float) (( *((unsigned char *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case SHORT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((short *) (fsrc + i * srcsize)) == SHORT_MAGIC) {
							*((float *) (fdst + i * dstsize)) = FLOAT_MAGIC;
						} else {
							if (bscaling) {
								*((float *) (fdst + i * dstsize)) = (float) (( *((short *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((float *) (fdst + i * dstsize)) = (float) (( *((short *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case INT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((int *) (fsrc + i * srcsize)) == INT_MAGIC) {
							*((float *) (fdst + i * dstsize)) = FLOAT_MAGIC;
						} else {
							if (bscaling) {
								*((float *) (fdst + i * dstsize)) = (float) (( *((int *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((float *) (fdst + i * dstsize)) = (float) (( *((int *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case FLOAT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((float *) (fsrc + i * srcsize)) == FLOAT_MAGIC) {
							*((float *) (fdst + i * dstsize)) = FLOAT_MAGIC;
						} else {
							if (bscaling) {
								*((float *) (fdst + i * dstsize)) = (float) (( *((float *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((float *) (fdst + i * dstsize)) = (float) (( *((float *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case DBL_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((double *) (fsrc + i * srcsize)) == DBL_MAGIC) {
							*((float *) (fdst + i * dstsize)) = FLOAT_MAGIC;
						} else {
							if (bscaling) {
								*((float *) (fdst + i * dstsize)) = (float) (( *((double *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((float *) (fdst + i * dstsize)) = (float) (( *((double *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				default:
					error_exit("convertarray: bad pixtype\n");
			}
			break;
		case DBL_PIXTYPE:
			dstsize = sizeof(double);
			switch(srcpixtype) {
				case UCHAR_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((unsigned char *) (fsrc + i * srcsize)) == UCHAR_MAGIC) {
							*((double *) (fdst + i * dstsize)) = DBL_MAGIC;
						} else {
							if (bscaling) {
								*((double *) (fdst + i * dstsize)) = (double) (( *((unsigned char *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((double *) (fdst + i * dstsize)) = (double) (( *((unsigned char *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case SHORT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((short *) (fsrc + i * srcsize)) == SHORT_MAGIC) {
							*((double *) (fdst + i * dstsize)) = DBL_MAGIC;
						} else {
							if (bscaling) {
								*((double *) (fdst + i * dstsize)) = (double) (( *((short *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((double *) (fdst + i * dstsize)) = (double) (( *((short *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case INT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((int *) (fsrc + i * srcsize)) == INT_MAGIC) {
							*((double *) (fdst + i * dstsize)) = DBL_MAGIC;
						} else {
							if (bscaling) {
								*((double *) (fdst + i * dstsize)) = (double) (( *((int *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((double *) (fdst + i * dstsize)) = (double) (( *((int *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case FLOAT_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((float *) (fsrc + i * srcsize)) == FLOAT_MAGIC) {
							*((double *) (fdst + i * dstsize)) = DBL_MAGIC;
						} else {
							if (bscaling) {
								*((double *) (fdst + i * dstsize)) = (double) (( *((float *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((double *) (fdst + i * dstsize)) = (double) (( *((float *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				case DBL_PIXTYPE:
					for (i = 0; i < nel; i++) {
						if (*((double *) (fsrc + i * srcsize)) == DBL_MAGIC) {
							*((double *) (fdst + i * dstsize)) = DBL_MAGIC;
						} else {
							if (bscaling) {
								*((double *) (fdst + i * dstsize)) = (double) (( *((double *) (fsrc + i * srcsize)) * bscale + bzero));
							} else {
								*((double *) (fdst + i * dstsize)) = (double) (( *((double *) (fsrc + i * srcsize))));
							}
						}
					}
					break;
				default:
					error_exit("convertarray: bad pixtype\n");
			}
			break;
		default:
			error_exit("convertarray: bad pixtype\n");
	}

}
