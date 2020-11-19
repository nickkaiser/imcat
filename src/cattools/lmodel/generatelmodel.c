/*
 * generatelmodel.c
 */


#define usage "\n\
NAME\n\
	generatelmodel --- generate lmodel for catalog of points\n\
\n\
SYNOPSIS\n\
	generatelmodel [-n aname] lmodel \n\
\n\
DESCRIPTION\n\
	'generatelmodel' reads from stdin a catalogue containing at least\n\
	a position vetor and a linear mode function superposition model\n\
	from 'lmodel' and computes the model function for each point in\n\
	the catalog.\n\
\n\
	The name and dimension of the position vector must match 'xname'\n\
	and 'xdim' in the lmodel header.\n\
\n\
	By default the realised values of the functions will be\n\
	named as the 'aname' header item in the lmodel header,\n\
	and will therefore overwrite any pre-existing opject item\n\
	of that name but you can supply an alternative with the -n option.\n\
\n\
	'generatelmodel' uses the imcat iostream mechanism, so you can supply\n\
	a paramater generating command 'pargencom |' in place of the lmodel argument.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "catlib/cat.h"
#include "utils/error.h"
#include "utils/args.h"
#include "utils/lmodel.h"
#include "utils/iostream.h"


main(int argc, char *argv[])
{
	/* catalogue stuff */
	cathead	*ipcat, *opcat, *lmodelcat;
	object	*ipobj, *opobj, *lmodelobj;
	item	*xitem, *aitem, *amoditem;
	int	xindex, aindex, amodindex, cattype;
	void	*xaddr, *amodaddr;

	/* model stuff */
	lmodel	*themodel;
	char	*amodname = NULL;
	double	fm, *x;
	void	*a;
	char	*lmodelname;
	iostream	*lmodelstream;
	int	m;

	/* args stuff */
	char	*flag;

	/* parse first 2 args */
	argsinit(argc, argv, usage);
	while (FLAG_ARG == nextargtype()) {
		flag = getflag();
		switch (flag[0]) {
			case 'n':
				amodname = getargs();
				break;
			default:
				error_exit(usage);
		}
	}
	lmodelname = getargs();

	/* read the lmodel */
	lmodelstream = openiostream(lmodelname, "r");
	if (!lmodelstream) {
		error_exit("generatelmodel: failed to open lmodel\n");
	}
	themodel = readlmodel(lmodelstream->f);

	/* set the name of the realised model */
	if (!amodname) {
		amodname = themodel->aname;
	}

	/* read input cat head */
	setcatipf(stdin);
	ipcat = readcathead();
	ipobj = newobject(ipcat);
        connectobjecttocathead(ipobj);
        allocobjectcontents(ipobj);

	/* get the xitem and check its type */
	xitem = getobjectitem(themodel->xname, ipcat);
	if (xitem->itype != NUM_TYPE || xitem->ndim != 1 || (xitem->dim)[0] != themodel->xdim) {
		error_exit("generatelmodel : mismatch between input cat x and lmodel x\n");
	}
	xindex = getobjectitemindex(themodel->xname, ipobj);

	/* create the op cat */
	opcat = (cathead *) calloc(1, sizeof(cathead));
	copyheaderinfo(opcat, ipcat);
	copycontentinfo(opcat, ipcat);
	amoditem = newitembydimarray(amodname, NUM_TYPE, (themodel->aitem)->ndim, (themodel->aitem)->dim);
	addobjectitem(amoditem, opcat);
	opobj = newobject(opcat);
	inheritcontents(opobj, ipobj);
	amodindex = getobjectitemindex(amodname, opobj);
	allocitemcontents(amoditem, &((opobj->addrlist)[amodindex]), 0);
	xaddr = (double *) ((ipobj->addrlist)[xindex]);
	amodaddr = (opobj->addrlist)[amodindex];

	/* add history string */
	addargscomment(argc, argv, opcat);

	/* set the type */
	getcatipfiletype(&cattype);
	setcatopfiletype(cattype);

	/* write header */
	writecathead(opcat);
	while (readobject(ipobj)) {
		zeromatrix(amodaddr, amoditem->ndim, amoditem->dim, 0);
		for (m = 0; m < themodel->nmodes; m++) {
			fm = lmodelfunc(themodel, m, xaddr);
			addtomatrix(amodaddr, (themodel->a)[m], fm, amoditem->ndim, amoditem->dim, 0);
		}
		writeobject(opobj);
	}
	exit(0);
}


