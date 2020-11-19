/*
 * getrange.c
 */

void	getrange(float *f, int npts, float *fmin, float *fmax)
{
	int	i;
	float	min, max;
	
	min = max = f[0];
	for (i = 1; i < npts; i++) {
		if (f[i] > max)
			max = f[i];
		if (f[i] < min)
			min = f[i];
	}
	*fmin = min;
	*fmax = max;
}
