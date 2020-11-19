 /*
 * pair.c
 */

#define usage "\n\n\
NAME\n\
	pair - generate a catalogue containing pairs of objects from input catalogue \n\
\n\
SYNOPSIS\n\
	pair [options...]\n\
		-d	# only output N (N - 1) / 2 distinct pairs\n\
\n\
DESCRIPTION\n\
        'pair' reads a catalogue from stdin and writes to stdout a catalogue\n\
	which has object items with the same names as in the input cat\n\
	but where each object is a 2-vector formed from a pair of\n\
	input objects. By default all N (N - 1) pairs are output.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n"


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "../utils/error.h"
#include "../catlib/cat.h"



main(int argc, char *argv[])
{
	int		arg = 1, distinct, i, dim[10], index;
	cathead		*ipcat, *opcat;
	object		*ipobj, *opobj, *ipobjbase = NULL, *obj1, *obj2;
	item		*ipitem, *opitem;
	void		**opaddr;

	/* defaults */
	distinct = 0;

	/* parse the args */
	while (arg < argc) {
		if (argv[arg][0] != '-')
			error_exit(usage);
		switch (argv[arg++][1]) {
			case 'd':
				distinct = 1;
				break;
			default:
				error_exit(usage);
				break;
		}
	}




	/* read the cat header */
	ipcat  = readcathead();
	

	/* create the output cathead */
	opcat = (cathead *) calloc(1, sizeof(cathead));
	copyheaderinfo(opcat, ipcat);

	/* add comment */
	addargscomment(argc, argv, opcat);


	/* create the object items */
	ipitem = ipcat->objectitembase;
	while (ipitem) {
		allocitemcontents(ipitem, &(ipitem->addr), 0);
		for (i = 0; i < ipitem->ndim; i++) {
			dim[i] = (ipitem->dim)[i];
		}
		opitem = newitem(ipitem->name, ipitem->itype, ipitem->ndim + 1, 2,
                        dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], dim[8], dim[9]);
                addobjectitem(opitem, opcat);
               	ipitem = ipitem->next;
	}

	/* and write output cathead */
	setcatopfiletype(BINARY_FILE_TYPE);
	writecathead(opcat);

	/* read the objects into a linked list */
	while (1) {
		ipobj = newobject(ipcat);
		allocobjectcontents(ipobj);
		if (!readobject(ipobj))
			break;
		ipobj->next = ipobjbase;
		ipobjbase = ipobj;
	}

	/* create the ouput object */
	opobj = newobject(opcat);
	allocobjectcontents(opobj);

	/* loop over pairs */
	obj1 = ipobjbase;
	while (obj1) {
		if (distinct) {
			obj2 = obj1;
		} else {
			obj2 = ipobjbase;
		}
		while (obj2) {
			if (obj2 != obj1) {
				for (index = 0; index < opobj->nitems; index++) {
					opaddr = (void **) (opobj->addrlist)[index];
					opaddr[0] = (obj1->addrlist)[index];
					opaddr[1] = (obj2->addrlist)[index];
				}
				writeobject(opobj);
			}
			obj2 = obj2->next;
		}
		obj1 = obj1->next;
	}
	exit(0);
}
