void	gaussian_kernel_filter(short **f, short **fs, int N1, int N2, int m, float rf);
void	kernel_filter(short **f, short **fs, int N1, int N2, int m, float (*filterfunc)(int di, int dj));
void	block_filter(short **f, short **fs, int N1, int N2, int m);
void	tukey(short **f, int N1, int N2);
void	schecterfilter(float **f, 
		int N1, int N2, 
		float **fs, 
		float sigma1, 
		float sigma2,
		float alpha,
		float magicsubstitute);
void	kolmogorovfilter(float **f, 
		int N1, int N2, 
		float **fs, 
		float sigma,
		float magicsubstitute);
void	gaussfilter(float **f, 
		int N1, int N2, 
		float **fs, 
		float A, 
		float B,
		float phi,
		float magicsubstitute);
void	mexicanfilter(float **f, 
		int N1, int N2, 
		float **fs, 
		float sigma1, 
		float sigma2,
		float magicsubstitute);
void	powerlawfilter(float **f, 
		int N1, int N2, 
		float **fs, 
		float alpha,
		float magicsubstitute);
void	acf(short **f, int N1, int N2, short **theacf, short magicsubstitute);
void	powerspectrum(float **f, int N1, int N2, float **P, int **nmodes);
float   gaussianfilterfunc(int i, int j);
float   schecterfilterfunction(float ki, float kj);
float   kolmogorovfilterfunction(float ki, float kj);
float   gaussballfunction(float ki, float kj);
float   gaussellipsoidfunction(float ki, float kj);
float	mexicanfilterfunction(float ki, float kj);
float	powerlawfilterfunction(float ki, float kj);


#ifndef PI
#define PI M_PI
#endif






