/*
 * gaussdev.c - from numerical recipes gasdev.c renamed to avod conflict with
 * various built in implementations 
 */

#include <stdlib.h>
#include <math.h>

double gaussdev(void)
{
        static int iset=0;
        static double gset;
        double fac,r,v1,v2;

        if  (iset == 0) {
                do {
                        v1=2.0*drand48()-1.0;
                        v2=2.0*drand48()-1.0;
                        r=v1*v1+v2*v2;
                } while (r >= 1.0);
                fac=sqrt(-2.0*log(r)/r);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}

void	seedgaussdev(int seed)
{
	srand48((long int) seed);
}
