/*
 * unmergecat.c
 */

#define usage "\n\n\
NAME\n\
	unmergecat - extract one member from a merged catalogue\n\
\n\
SYNOPSIS\n\
	unmergecat N\n\
\n\
DESCRIPTION\n\
	'unmergecat' reads a merged catalogue created by 'mergecats'.\n\
	Ignorring any object items which are not vectors of things\n\
	(e.g. the 'detection mask' or number of detections column)\n\
	it creates a new catalogue in which each object item\n\
	is the (N-1)th element of the item of the same name in the\n\
	input catalogue.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n"


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "../../utils/error.h"
#include "../../catlib/cat.h"


main(int argc, char *argv[])
{
	int		arg, N, first, firstdim0, i, dim[10];
	cathead		*ipcat, *opcat;
	object		*ipobj, *opobj;
	item		*ipitem, *opitem;

	/* handle args */
	if ((argc != 2) || ((argc ==2) && (argv[1][0] == '-')))
		error_exit(usage);
	if (1 != sscanf(argv[1], "%d", &N))
		error_exit(usage);

	/* read the input cat headers */
	ipcat  = readcathead();

	/* create the output cathead */
	opcat = (cathead *) calloc(1, sizeof(cathead));
	copyheaderinfo(opcat, ipcat);

	/* add comment, headeritems */
	addargscomment(argc, argv, opcat);

	/* copy appropriate element of vector type objects*/
	ipitem = ipcat->objectitembase;
	first = 1;
	while (ipitem) {
		if (ipitem->ndim == 1) {	/* skip the maks, ndet columns if present */
			allocitemcontents(ipitem, &(ipitem->addr), 0);
			ipitem = ipitem->next;
			continue;
		}
		if (first) {
			firstdim0 = (ipitem->dim)[0];
			if (firstdim0 < N)
				error_exit("unmergecat: less than N elements in object item\n");
			first = 0;
		}
		if ((ipitem->dim)[0] != firstdim0)
			error_exit("unmergecat: object items do not have identical dim[0]!\n");
		for (i = 0; i < ipitem->ndim - 1; i++) {
			dim[i] = (ipitem->dim)[i + 1];
		}
		opitem = newitem(ipitem->name, ipitem->itype, ipitem->ndim - 1,
			dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], dim[8], dim[9]);
		addobjectitem(opitem, opcat);
		allocitemcontents(ipitem, &(ipitem->addr), 0);
		opitem->addr = ((void **)(ipitem->addr))[N];
		ipitem = ipitem->next;
	}

	/* create the objects */
	ipobj = newobject(ipcat);
	connectobjecttocathead(ipobj);	
	opobj = newobject(opcat);
	connectobjecttocathead(opobj);	

	/* and write it out */
	setcatopfiletype(BINARY_FILE_TYPE);
	writecathead(opcat);
	while (readobject(ipobj)) {
		writeobject(opobj);
	}
	exit(0);
}



