/*
 * getsky.h
 */



int	getsky(double *fb0, double *dfb, double r1, double r2, double *x, 
		float **f, int N1, int N2);

void	pruneextreme(void);
int	getplane(double *fb0, double *dfb);
int	fmodecmp(const void *ptr1, const void *ptr2);
