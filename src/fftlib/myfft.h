/*
 * myfft.h
 *
 *
 * general comments:
 *
 *	the problem is that different FFT routines pack their transforms
 *	in different ways, and some only work for image dimensions 2^N
 *
 *	solution here is to implement a set of functions to
 *	do filtering, power spectrum analysis, auto/cross-correlation
 *	analysis in such a way that the calling routine only
 *	specifies the real image size, never has to figure out the
 *	type of transform data, its size, or where
 *	a particular frequency lives in the transform
 *
 *	2-D (float **) images f[y][x] only
 *	fast dimension Nx (Nx=N1 in fits notation)
 *
 *
 * this file declares the following functions (which are implemented in
 * the various platform secific fft_xxx.c files):
 * 	forward_fft(f, Nx, Ny, fk)
 *	inverse_fft(fk, Nx, Ny, f)
 *	filter(fk, Nx, Ny, filterfunc)
 *	cfilter(fk, Nx, Ny, rfunc, ifunc)
 * 	ccf(fk1, fk2, Nx, Ny, ff, x0, y0)
 * 	power(fk, Nx, Ny, p, x0, y0)
 * 	alloc_fft(&fk, Nx, Ny)
 * 	free_fft(fk, Nx, Ny)
 * 	get_fft(fk, Nx, Ny, fdest)
 * 	set_fft(fk, Nx, Ny, fsrc)
 *
 *
 * forward_fft(f, Nx, Ny, &fk)
 *	takes a imtype array f[y][x] and creates a (float *)
 *	array fk containing (half) the complex transform.
 *
 *
 * inverse_fft(fk, Nx, Ny, f)
 *	does inverse transform
 *
 *
 * filter(fk, Nx, Ny, filterfunc)
 *	applies a smoothing filter filterfunc(kx, ky)
 * 
 *
 * cfilter(fk, Nx, Ny, rfunc, ifunc)
 *	applies a complex filter f_k => f_k * (rfunc + i * ifunc)
 * 
 *
 * ccf(fk1, fk2, Nx, Ny, ff, x0, y0)
 *	calculates cross-correlation function ff12[y][x] of two
 *	Nx x Ny images with tranforms fk1[y][x], fk2[y][x]
 *	with origin shifted to x0, y0
 * 
 *
 * power(fk, Nx, Ny, p, x0, y0)
 *	calculates power spectrum p[ky][[kx] from transform fk
 *	with origin shifted to kx0, ky0
 *	
 *	
 * alloc_fft(fk, Nx, Ny)
 *	allocates memory for transform
 *	
 *	
 * alloc_fft_shm(fk, Nx, Ny)
 *	allocates memory for transform - shared memory version
 *	
 *	
 * free_fft(fk, Nx, Ny)
 *	frees memory for transform
 *	
 *	
 * free_fft_shm(fk, Nx, Ny, shmid)
 *	frees memory for transform - shared memory version
 *	
 *	
 * copy_fft(fksrc, Nx, Ny, fkdest)
 *	makes a copy fkdest of fksrc
 *	
 *	
 * get_fft(fk, Nx, Ny, fdest)
 *	fill the 2-d float array fdest with the real and imaginary
 *	parts of fk.  The (y * x) dimensions of fdest are (2 Ny) * Nx.
 *	The real parts are stored in the Ny * Nx array with
 *	base fdest and with zero-frequency at Ny / 2, Nx / 2.
 *	The imaginary parts are stored in the Ny * Nx array
 *	with base (fdest + Ny).  Only works if Nfftx = Nx, Nffty = Ny.
 *
 * set_fft(fk, Nx, Ny, fsrc)
 *	copy fk values from (2 * Ny) * Nx array with storage
 *	convention as for get_fft().  Only works if Nfftx = Nx, Nffty = Ny.
 *
 * The FFTPACK, FFTW versions also define complex fft's
 * the format here is that fsrc, fdst are size (2 ny) * nx
 * and contain concatenation of real and imaginary images.
 * zero frequency still lives at nx / 2, ny / 2 in each sub-image
 *	forward_cfft(float **fsrc, int nx, int ny, float **fdst)
 *	inverse_cfft(float **fsrc, int nx, int ny, float **fdst)
 *	
 * In addition there are some functions defined in my fft.c	
 *	
 * cycleimage(f, Nx, Ny, x0, y0)
 *	shifts the origin of a periodic image from 0,0 to x0,y0	
 *	
 *	
 * substitute(f, Nx, Ny, magicsubstitute)
 *	replaces MAGIC values
 *
 */



#ifdef SGI_FFT
#include "fft_SGI.h"
#endif

#ifdef DEC_FFT
#include "fft_DEC.h"
#endif

#ifdef FFTPACK_FFT
#include "fft_FFTPACK.h"
#endif

#ifdef NR_FFT
#include "fft_NR.h"
#endif

#ifdef FFTW_FFT
#include <fftw.h>
#include "fft_FFTW.h"
#endif


#define SHORT_IMAGE	0
#define FLOAT_IMAGE	1
#define DOUBLE_IMAGE	2


void	forward_fft(float **f, int Nx, int Ny, fft_type fk);
void	inverse_fft(fft_type fk, int Nx, int Ny, float **f);
void	filter(fft_type fk, int Nx, int Ny, float (*filterfunc)(float ki, float kj));
void	cfilter(fft_type fk, int Nx, int Ny, float (*rfunc)(float ki, float kj),  float (*ifunc)(float ki, float kj));
void	ccf(fft_type fk1, fft_type fk2, int Nx, int Ny, float **ff12, int x0, int y0);
void	power(fft_type fk, int Nx, int Ny, float **p, int x0, int y0);
void	cycleimage(float **f, int Nx, int Ny, int x0, int y0);
void	alloc_fft(fft_type *fk, int Nx, int Ny);
void	free_fft(fft_type fk, int Nx, int Ny);
#ifdef DUAL_PROC
int	alloc_fft_shm(fft_type *fk, int Nx, int Ny);
void	free_fft_shm(fft_type fk, int shmid);
void	set_shmids(int fr_shmid, int fk_shmid);
void	set_frshmid(int fr_shmid);
void	set_fkshmid(int fk_shmid);
void	get_shmids(int *fr_shmid, int *fk_shmid);
#endif
void	substitute(float **f, int Nx, int Ny, float magicsubstitute);
void	copy_fft(fft_type fksrc, int Nx, int Ny, fft_type fkdest);
void	get_fft(fft_type fk, int Nx, int Ny, float **fdest);
void	set_fft(fft_type fk, int Nx, int Ny, float **fsrc);
 
void	forward_cfft(float **fsrc, int nx, int ny, float **fdst);
void	inverse_cfft(float **fsrc, int nx, int ny, float **fdst);
