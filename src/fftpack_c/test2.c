#include <stdio.h>

main(int argc, char *argv[])
{
	int	n;
	float	*r, *work;

	cffti_(&n, work);
	cfftf_(&n, r, work);
	cfftb_(&n, r, work);

	rffti_(&n, work);
	rfftf_(&n, r, work);
	rfftb_(&n, r, work);
}
