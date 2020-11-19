/*
 * fft_FFTPACK.h
 */

typedef struct { float r, i; } complex;

int     cffti_(int *n, float *work);
int     cfftf_(int *n, complex *c, float *work);
int     cfftb_(int *n, complex *c, float *work);

int     rffti_(int *n, float *work);
int     rfftf_(int *n, float *r, float *work);
int     rfftb_(int *n, float *r, float *work);

int     ezffti_(int *n, float *work);
int     ezfftf_(int *n, float *r, float *a0, float *a, float *b, float *work);
int     ezfftb_(int *n, float *r, float *a0, float *a, float *b, float *work);

typedef		complex **fft_type;
