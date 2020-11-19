/*
 * fitstack_modefunc.c
 */


double	f(int l, int m, double *x)
{
	double	res = 1;
	int	i;

	for (i = 0; i < (l - m) ; i++) {
		res *= x[0];
	}
	for (i = 0; i < m ; i++) {
		res *= x[1];
	}
	return (res);
}
