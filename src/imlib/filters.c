/**
 **  filters.c --- image processing routines
 **/
 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>

#include "../utils/error.h"
#include "filters.h"
#include "fits.h"
#include "../fftlib/myfft.h"

#define	max(a, b) ((a) < (b) ? (b) : (a))
#define	min(a, b) ((a) < (b) ? (a) : (b))

/* globals for the various kernel functions */
static	float	gaussballsigma;
static	float	gaussfilterparam;
static	float	schectersigma1, schectersigma2, schecteralpha;
static	float	kolmogorovsigmasquared;
static	float	exponentialsigmasquared, exponentialgamma;
static	float	gausssigma11, gausssigma22, gausssigma12;
static	float	mexicansigma1, mexicansigma2;
static	float	powerlawalpha;

void	gaussian_kernel_filter(float **f, float **fs, int N1, int N2, int m, float rf)
{
	gaussfilterparam = 0.5 / (rf * rf);
	kernel_filter(f, fs, N1, N2, m, gaussianfilterfunc);
}


void	kernel_filter(float **f, float **fs, int N1, int N2, int m, float (*filterfunc)(int di, int dj))
{
	int 		i, j, ii, jj, bad;
	long		result;
	float 		**k, k_sum = 0.0, fsum;
	
	if (2 * (m / 2) == m)					/* create kernel */
		error_exit("kernel_filter: kernel size must be odd");
	k = (float **) calloc(m, sizeof(float));
	if (!k)
		error_exit("kernel_filter: memory allocation failure");
	k += m / 2;
	for (i = - m / 2; i <= m / 2; i++) {
		k[i] = (float *) calloc(m, sizeof(float));
		if (!k[i])
			error_exit("kernel_filter: memory allocation failure");
		k[i] += m / 2;
		for (j = - m / 2; j <= m / 2; j++) {
			k[i][j] = filterfunc(i, j);
			k_sum += k[i][j];
		}
	}
			
	for (i = 0; i < N2; i++) {				/* do the filtering */
		for (j = 0; j < N1; j++) {
			bad = result = 0;
			fsum = 0.0;
			for (ii = - m / 2; ii <= m / 2; ii++) {
				if (i + ii >= N2 || i + ii < 0)
					continue;
				for (jj = - m / 2; jj <= m / 2; jj++) {
					if (j + jj >= N1 || j + jj < 0)
						continue;
					if (f[i+ii][j+jj] == SHORT_MAGIC)
						bad = 1;
					fsum += k[ii][jj] * f[i + ii][j + jj];
				}
			}
			result = floor(0.5 + fsum /  k_sum);
			fs[i][j] = (float) max(SHRT_MIN, min(SHRT_MAX, result));			
			if (bad)
				fs[i][j] = SHORT_MAGIC;
		}
	}
	
	for (i = - m / 2; i <= m / 2; i++) 		/* destroy kernel */
		free(k[i] - m / 2);
	free(k - m / 2);
}


void	block_filter(float **f, float **fs, int N1, int N2, int m)
{
	int 		i, j, ii, jj, bad;
	float		result;
	
	if (2 * (m / 2) == m)
		error_exit("block_filter: block size must be odd");
			
	for (i = 0; i < N2; i++) {				/* do the filtering */
		for (j = 0; j < N1; j++) {
			bad = 0;
			result = 0.0;
			for (ii = - m / 2; ii <= m / 2; ii++) {
				if (i + ii >= N2 || i + ii < 0)
					continue;
				for (jj = - m / 2; jj <= m / 2; jj++) {
					if (j + jj >= N1 || j + jj < 0)
						continue;
					if (f[i+ii][j+jj] == FLOAT_MAGIC)
						bad = 1;
					result += f[i + ii][j + jj];
				}
			}
			if (bad) {
				fs[i][j] = FLOAT_MAGIC;
			} else {
				fs[i][j] = result / (m * m);
			}
		}
	}
}





void		tukey(float **f, int N1, int N2)
{
	int 	i, j;
	float	f_last;
	
	for (i = 1; i < N2 - 1; i++) {
		f_last = f[i][0];
		for (j = 1; j < N1 - 1; j++) {
			if (f[i][j] > f_last && f[i][j] > f[i][j+1])
				f[i][j] = max(f_last, f[i][j+1]);
			if (f[i][j] < f_last && f[i][j] < f[i][j+1])
				f[i][j] = min(f_last, f[i][j+1]);
			f_last = f[i][j];
		}
	}
	for (j = 1; j < N1 - 1; j++) {
		f_last = f[0][j];
		for (i = 1; i < N2 - 1; i++) {
			if (f[i][j] > f_last && f[i][j] > f[i+1][j])
				f[i][j] = max(f_last, f[i+1][j]);
			if (f[i][j] < f_last && f[i][j] < f[i+1][j])
				f[i][j] = min(f_last, f[i+1][j]);
			f_last = f[i][j];
		}
	}
}



/*
 * implements (1 + k^2 sigma1^1)^(-alpha/2) exp(-0.5 k^2 sigma2^2)
 * smoothing
 */
void	schecterfilter(float **f, 
		int N1, int N2,
		float **fs, 
		float sigma1, 
		float sigma2,
		float alpha,
		float magicsubstitute)
{
	fft_type 	fk;
	int		i, j;
#ifdef DUAL_PROC
	int		fk_shmid;
#endif

	schectersigma1 = sigma1;
	schectersigma2 = sigma2;
	schecteralpha = alpha;

	/* if fs != f we want to make a copy to fs and transform that instead */
	/* so that f is preserved intact */
	if (fs != f) {
		for (i = 0; i < N2; i++) {
			for (j = 0; j < N1; j++) {
				fs[i][j] = f[i][j];
			}
		}
	}
	substitute(fs, N1, N2, magicsubstitute);
#ifdef DUAL_PROC
	fk_shmid = alloc_fft_shm(&fk, N1, N2);
	set_fkshmid(fk_shmid); 
#else
	alloc_fft(&fk, N1, N2);
#endif
	forward_fft(fs, N1, N2, fk);
	filter(fk, N1, N2, schecterfilterfunction); 
	inverse_fft(fk, N1, N2, fs);
#ifdef DUAL_PROC
	free_fft_shm(fk, fk_shmid);
#else
	free_fft(fk, N1, N2);
#endif
}




/*
 * implements exp(-0.5 (k sigma)^(5/3))
 * smoothing
 */
void	kolmogorovfilter(float **f, 
		int N1, int N2,
		float **fs, 
		float sigma,
		float magicsubstitute)
{
	fft_type 	fk;
	int		i, j;
#ifdef DUAL_PROC
	int		fk_shmid;
#endif

	kolmogorovsigmasquared = sigma * sigma;

	/* if fs != f we want to make a copy to fs and transform that instead */
	/* so that f is preserved intact */
	if (fs != f) {
		for (i = 0; i < N2; i++) {
			for (j = 0; j < N1; j++) {
				fs[i][j] = f[i][j];
			}
		}
	}
	substitute(fs, N1, N2, magicsubstitute);
#ifdef DUAL_PROC
	fk_shmid = alloc_fft_shm(&fk, N1, N2);
	set_fkshmid(fk_shmid); 
#else
	alloc_fft(&fk, N1, N2);
#endif
	forward_fft(fs, N1, N2, fk);
	filter(fk, N1, N2, kolmogorovfilterfunction); 
	inverse_fft(fk, N1, N2, fs);
#ifdef DUAL_PROC
	free_fft_shm(fk, fk_shmid);
#else
	free_fft(fk, N1, N2);
#endif
}

/*
 * implements exp(-(k sigma)^(gamma))
 * smoothing
 */
void	exponentialfilter(float **f, 
		int N1, int N2,
		float **fs, 
		float sigma,
		float gamma,
		float magicsubstitute)
{
	fft_type 	fk;
	int		i, j;
#ifdef DUAL_PROC
	int		fk_shmid;
#endif

	exponentialsigmasquared = sigma * sigma;
	exponentialgamma = gamma;

	/* if fs != f we want to make a copy to fs and transform that instead */
	/* so that f is preserved intact */
	if (fs != f) {
		for (i = 0; i < N2; i++) {
			for (j = 0; j < N1; j++) {
				fs[i][j] = f[i][j];
			}
		}
	}
	substitute(fs, N1, N2, magicsubstitute);
#ifdef DUAL_PROC
	fk_shmid = alloc_fft_shm(&fk, N1, N2);
	set_fkshmid(fk_shmid); 
#else
	alloc_fft(&fk, N1, N2);
#endif
	forward_fft(fs, N1, N2, fk);
	filter(fk, N1, N2, exponentialfilterfunction); 
	inverse_fft(fk, N1, N2, fs);
#ifdef DUAL_PROC
	free_fft_shm(fk, fk_shmid);
#else
	free_fft(fk, N1, N2);
#endif
}


void	gaussfilter(float **f, 
		int N1, int N2, 
		float **fs, 
		float A, 
		float B,
		float phi,
		float magicsubstitute)
{
	fft_type 	fk;
	float	c, s;
	int		i, j;
#ifdef DUAL_PROC
	int		fk_shmid;
#endif
	
	if (A == B) {
		gaussballsigma = A * A;
	} else {
		phi *= 2 * PI / 360;
		c = cos(phi);
		s = sin(-phi);
		gausssigma11 = A * A * c * c + B * B * s * s;
		gausssigma22 = A * A * s * s + B * B * c * c;
		gausssigma12 = (B * B - A * A) * c * s;
	}

	/* if fs != f we want to make a copy to fs and transform that instead */
	/* so that f is preserved intact */
	if (fs != f) {
		for (i = 0; i < N2; i++) {
			for (j = 0; j < N1; j++) {
				fs[i][j] = f[i][j];
			}
		}
	}
	substitute(fs, N1, N2, magicsubstitute);
#ifdef DUAL_PROC
	fk_shmid = alloc_fft_shm(&fk, N1, N2);
	set_fkshmid(fk_shmid); 
#else
	alloc_fft(&fk, N1, N2);
#endif
	forward_fft(fs, N1, N2, fk);
	if (A == B) {
		filter(fk, N1, N2, gaussballfunction); 
	} else {
		filter(fk, N1, N2, gaussellipsoidfunction); 
	}
	inverse_fft(fk, N1, N2, fs);
#ifdef DUAL_PROC
	free_fft_shm(fk, fk_shmid);
#else
	free_fft(fk, N1, N2);
#endif
}



void	mexicanfilter(float **f, 
					int N1, int N2,
					float **fs, 
					float sigma1, 
					float sigma2,
					float magicsubstitute)
{
	fft_type 	fk;
#ifdef DUAL_PROC
	int		fk_shmid;
#endif

	mexicansigma1 = sigma1;
	mexicansigma2 = sigma2;

	substitute(f, N1, N2, magicsubstitute);
#ifdef DUAL_PROC
	fk_shmid = alloc_fft_shm(&fk, N1, N2);
	set_fkshmid(fk_shmid); 
#else
	alloc_fft(&fk, N1, N2);
#endif
	forward_fft(f, N1, N2, fk);
	filter(fk, N1, N2, mexicanfilterfunction); 
	inverse_fft(fk, N1, N2, fs);
#ifdef DUAL_PROC
	free_fft_shm(fk, fk_shmid);
#else
	free_fft(fk, N1, N2);
#endif

}


void	powerlawfilter(float **f, 
					int N1, int N2,
					float **fs, 
					float alpha,
					float magicsubstitute)
{
	fft_type 	fk;
#ifdef DUAL_PROC
	int		fk_shmid;
#endif

	powerlawalpha = alpha;

	substitute(f, N1, N2, magicsubstitute);
#ifdef DUAL_PROC
	fk_shmid = alloc_fft_shm(&fk, N1, N2);
	set_fkshmid(fk_shmid); 
#else
	alloc_fft(&fk, N1, N2);
#endif
	forward_fft(f, N1, N2, fk);
	filter(fk, N1, N2, powerlawfilterfunction); 
	inverse_fft(fk, N1, N2, fs);
#ifdef DUAL_PROC
	free_fft_shm(fk, fk_shmid);
#else
	free_fft(fk, N1, N2);
#endif

}


void	powerspectrum(float **f, int N1, int N2, float **P, int **nmodes)
/** 2-D power spectrum analysis. f[][] is N2 * N1 array. N must be power of two.
 ** P[] and nmodes[] must be allocated by calling function with dimension N / 2.
 **
 **	Result is P[k] = sum f exp(ik.r) averages on rings in k-space.
 ** Array nmodes[k] contains the number of modes used in averaging which
 ** is needed for uncertainty.
 **/
{
	error_exit("powerspectrum: to be implemented!\n");
}


#define ZMAX	10

float	gaussianfilterfunc(int i, int j)
{
	float	z;

	z = gaussfilterparam * (i * i + j * j);
        if (z > ZMAX)
		return (0.0);
	else
		return (exp(- gaussfilterparam * (i * i + j * j)));
}


float	schecterfilterfunction(float ki, float kj)
{
	float	kk;
	
	kk = ki * ki + kj * kj;
	return(pow(1 + kk * schectersigma1 * schectersigma1, - 0.5 * schecteralpha) *
				exp(-0.5 * kk * schectersigma2 * schectersigma2));
}


float	kolmogorovfilterfunction(float ki, float kj)
{
	float	kk;
	
	kk = ki * ki + kj * kj;
	return(exp(-0.5 * pow(kk * kolmogorovsigmasquared, 0.833333)));
}

float	exponentialfilterfunction(float ki, float kj)
{
	float	kk;
	
	kk = ki * ki + kj * kj;
	return(exp(- pow(kk * exponentialsigmasquared, 0.5 * exponentialgamma)));
}


float	gaussballfunction(float ki, float kj)
{
	float	z;

	z = 0.5 * gaussballsigma * (ki * ki +  kj * kj);
        if (z > ZMAX)
		return (0.0);
	else
		return (exp(- z));
}



float	gaussellipsoidfunction(float ki, float kj)
{
	float	z;

	z = 0.5 * (gausssigma11 * ki * ki + 2 * gausssigma12 * ki * kj +
		gausssigma22 * kj * kj);
        if (z > ZMAX)
		return (0.0);
	else
		return (exp(- z));
}



float	mexicanfilterfunction(float ki, float kj)
{
	float	kk;
		
	kk = ki * ki + kj * kj;
	return (exp(-0.5 * kk * mexicansigma1 * mexicansigma1) -
		exp(-0.5 * kk * mexicansigma2 * mexicansigma2));
}


float	powerlawfilterfunction(float ki, float kj)
{
	float	kk;
		
	kk = ki * ki + kj * kj;
	return (kk != 0.0 ? pow(kk, powerlawalpha / 2) : 0.0);
}



#undef	ZMAX

