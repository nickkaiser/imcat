/*
 * apphot.h
 */

typedef struct pixel {
	double	r;
	double	f;
	double	fcum;
	double	nu;
	double	nucum;
	double	nus;
} pixel;

int	apphot(double *flux, double *rh, double *rql, double *rqu, double *rp, 
		double *nbad, double *fmax,
		double *x, double rap, double *fb0, double *dfb, 
		float **f, int N1, int N2);


double	rpetrosian(double *x, float **f, int N1, int N2, double *fb0, double *dfb);

int	pixcmp(pixel *pix1, pixel *pix2);
