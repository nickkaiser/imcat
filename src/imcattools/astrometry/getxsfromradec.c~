/*
 * getxsfromradec.c
 */

#define usage "\n\
NAME\n\
	getxsfromradec --- convert from ra, dec to stereographic x-coord\n\
\n\
SYNOPSIS\n\
	getxsfromradec [-i | -u] RA0 DEC0 [THETA]\n\
\n\
DESCRIPTION\n\
	getxsfromradec reads from stdin a catalogue containing celestial\n\
	coordinates (and perhaps other stuff) and computes the\n\
	stereographic projection xs with tangent point \n\
	The orientation is such that xs[0] increases with decreasing\n\
	RA, and x[1] increases with DEC.\n\
	The resulting vector is expressed in degrees.\n\
\n\
	Tangent point args may be decimal (in degrees) or colon separated\n\
	triplets interpreted as h:m:s for ra and d:m:s for dec.\n\
\n\
	If the third optional argument THETA is given then the\n\
	stereographic coords will be rotated anticlockwise through\n\
	THETA degrees. THETA must be given in decimal notation.\n\
\n\
	With the -i flag we perform the inverse operation.\n\
\n\
AUTHOR\n\
	Nick Kaiser	kaiser@hawaii.edu\n\
\n\n"


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


#include "catlib/cat.h"
#include "utils/error.h"
#include "utils/args.h"
#include "getxfromradecfunc.h"
#include "radecio.h"

main (int argc, char *argv[])
{
	double	*RA, *DEC, RA0, DEC0, *x, theta, c, s, xp[2];
	cathead	*inputcathead, *outputcathead;
	item	*xitem, *RAitem, *DECitem;
	int	ix, iRA, iDEC, angletype, filetype, rotate, doinverse;
	object	*inputobject, *outputobject;
	char	*flag;	

	/* defaults */
	doinverse = 0;

	argsinit(argc, argv, usage);
	while (nextargtype() == FLAG_ARG) {
		flag = getflag();
		switch (flag[0]) {
			case 'i':
				doinverse = 1;
				break;
			default:
				error_exit(usage);
				break;
		}
	}
	RA0  = getra(getargs(), &angletype);
	DEC0 = getdec(getargs(), &angletype);
	if (nextargtype()) {
		rotate = 1;
		theta = getargd();
		c = cos(M_PI * theta / 180.0);
		s = sin(M_PI * theta / 180.0);
	} else {
		rotate = 0;
	}

	inputcathead = readcathead();
	inputobject = newobject(inputcathead);
	connectobjecttocathead(inputobject);
	allocobjectcontents(inputobject);
	if (doinverse) {
		x  =  (double *) ((inputobject->addrlist)[getobjectitemindex("xs",  inputobject)]);
	} else {
		RA  = (double *) ((inputobject->addrlist)[getobjectitemindex("RA",  inputobject)]);
		DEC = (double *) ((inputobject->addrlist)[getobjectitemindex("DEC", inputobject)]);
	}

	outputcathead = (cathead *) calloc(1, sizeof(cathead));	
	copyheaderinfo(outputcathead, inputcathead);
	addargscomment(argc, argv, outputcathead);
	copycontentinfo(outputcathead, inputcathead);
	if (doinverse) {
		RAitem = newitem("RA", NUM_TYPE, 1, 1);
		DECitem = newitem("DEC", NUM_TYPE, 1, 1);
		addobjectitem(RAitem, outputcathead);
		addobjectitem(DECitem, outputcathead);
		outputobject = newobject(outputcathead);
		iRA = getobjectitemindex("RA", outputobject);
		iDEC = getobjectitemindex("DEC", outputobject);
		inheritcontents(outputobject, inputobject);
		allocitemcontents(RAitem, &((outputobject->addrlist)[iRA]), 0);
		allocitemcontents(DECitem, &((outputobject->addrlist)[iDEC]), 0);
		RA  =  (double *) ((outputobject->addrlist)[iRA]);
		DEC  =  (double *) ((outputobject->addrlist)[iDEC]);
	} else {
		xitem = newitem("xs", NUM_TYPE, 1, 2);
		addobjectitem(xitem, outputcathead);
		outputobject = newobject(outputcathead);
		ix = getobjectitemindex("xs", outputobject);
		inheritcontents(outputobject, inputobject);
		allocitemcontents(xitem, &((outputobject->addrlist)[ix]), 0);
		x  =  (double *) ((outputobject->addrlist)[ix]);
	}
	getcatipfiletype(&filetype);
	setcatopfiletype(filetype);
	writecathead(outputcathead);

	while (readobject(inputobject)) { 
		if (doinverse) {
			if (rotate) {
				xp[0] =  c * x[0] + s * x[1];
				xp[1] =  -s * x[0] + c * x[1];
			} else {
				xp[0] = x[0];
				xp[1] = x[1];
			}
			inversegetxcoords(RA, DEC, RA0, DEC0, xp);
		} else {
			getxcoords(*RA, *DEC, RA0, DEC0, x);
			if (rotate) {
				xp[0] = x[0];
				xp[1] = x[1];
				x[0] =  c * xp[0] - s * xp[1];
				x[1] =  s * xp[0] + c * xp[1];
			}
		}
		writeobject(outputobject);
	}

	exit(0);
}
