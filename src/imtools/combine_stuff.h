/*
 * combine_stuff.h
 */


float   	rankedimage(float *f, int nel, int rank);
float   	median(float *f, int nel);
float   	getsigma(float *sigmavec, int nel);
float   	mean(float *f, int nel);
float   	avsigclip(float *f, int nel, float clip);
float   	avsigclip2(float *f, float *sigma, int nel, float clip, float frac1 , float frac2,
			float *sigmaout);
