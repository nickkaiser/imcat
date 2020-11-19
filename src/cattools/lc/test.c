#include <stdio.h>

#include "../../catlib/cat.h"

main(int argc, char *argv[]){
	cathead	*inputcathead, *outputcathead;
	object	*inputobject, *outputobject;
	item	*theitem, *textitem;
	int	xindex, findex, filenameindex, imindex;
	float	*x;
	double	**g, *d, *f;
	char	**filename, text[128], **im;

	inputcathead = readcathead();

	g = (double **) getheaderitemaddress("G", inputcathead);
	g[0][0] = -100.0;
	d = (double *) calloc(2, sizeof(double));
	d[0] = 1.3;
	d[1] = 5.7;
	setheaderitemaddress("D", inputcathead, (void *) d);


	outputcathead = copycathead(inputcathead);
 	addcomment("F comes from d", outputcathead);
	addargscomment(argc, argv, outputcathead);
	deleteheaderitem("WORDS", outputcathead);
	theitem = newitem("F", NUM_TYPE, 1, 2);
	addobjectitem(theitem, outputcathead);
	deleteobjectitem("IM", outputcathead);
	textitem = newitem("FILENAME", TEXT_TYPE, 1, 1);
	addobjectitem(textitem, outputcathead);
	writecathead(outputcathead);


	inputobject = newobject(inputcathead);
	allocobjectcontents(inputobject);

	xindex = getobjectitemindex("X", inputobject);
	x = (float *) getaddress(inputobject, xindex);

	imindex = getobjectitemindex("IM", inputobject);
	im = (char **) getaddress(inputobject, imindex);


	outputobject = newobject(outputcathead);
	inheritcontents(outputobject, inputobject);
	allocobjectcontents(outputobject);

	findex = getobjectitemindex("F", outputobject);
/*	f = (float *) getaddress(outputobject, findex);*/
	setaddress(outputobject, findex, (void *) d);

	filenameindex = getobjectitemindex("FILENAME", outputobject);
	filename = (char **) getaddress(outputobject, filenameindex);

	while (readobject(inputobject)) {
		strcpy(text, im[0]);
		strcat(text, " and ");
		strcat(text, im[1]);
		copystring(&(filename[0]), text);
/*		f[0] = 2.0 * x[0];
		f[1] = 3.0 * x[0];*/
		writeobject(outputobject);
	}
}
