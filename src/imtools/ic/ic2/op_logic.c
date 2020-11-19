/*
 * op_logic.c
 */


double	gt(double x, double y)
{
	return ((double) (x > y));
}



double	ge(double x, double y)
{
	return ((double) (x >= y));
}



double	lt(double x, double y)
{
	return ((double) (x < y));
}



double	le(double x, double y)
{
	return ((double) (x <= y));
}



double	ne(double x, double y)
{
	return ((double) (x != y));
}


double	not(double x)
{
	return ((double) (x == 0.0));
}

