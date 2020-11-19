/*
 * op_rand.c - random number generator
 */

#include <stdio.h>
#include "../../catlib/cat.h"
#include "stack.h"
#include "error.h"
#include "math.h"

double	drand48(void);
double gasdev(void);

void	randinit(item *operator)
{
	item	*result;
	double	*randno;

	result = operator->next = newitem("RAND_TEMP", NUM_TYPE, 1, 1);
	randno = (double *) calloc(1, sizeof(double));
	result->addr = (void *) randno;
	*randno = drand48();
}

void	randdoit(item *operator)
{
	item	*result;
	double	*randno;

	randno = (double *) ((operator->next)->addr);
	*randno = drand48();
}

void	grandinit(item *operator)
{
	item	*result;
	double	*randno;

	result = operator->next = newitem("RAND_TEMP", NUM_TYPE, 1, 1);
	randno = (double *) calloc(1, sizeof(double));
	result->addr = (void *) randno;
	*randno = gasdev();
}

void	granddoit(item *operator)
{
	item	*result;
	double	*randno;

	randno = (double *) ((operator->next)->addr);
	*randno = gasdev();
}


double gasdev(void)
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
