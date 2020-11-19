#include <stdio.h>
#include "catlib/cat.h"
#include "lmodel.h"

main (int argc, char *argv[]) 
{
	lmodel	*srcmodel, *dstmodel;

	srcmodel = readlmodel(stdin);


	writelmodel(srcmodel, stdout);
/**/
/*

	dstmodel = newfourierlmodel(srcmodel->xitem, srcmodel->aitem, 0, 3, 1, 1.0);
	dstmodel = newzernikelmodel(srcmodel->xitem, srcmodel->aitem, 0, 3);

	writelmodel(dstmodel, stdout);
*/
	

}
