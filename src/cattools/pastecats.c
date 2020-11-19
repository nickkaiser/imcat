/*
 * pastecats.c
 */

#define usage "\n\n\
NAME\n\
	pastecats - combine catalogues \n\
\n\
SYNOPSIS\n\
	pastecats a.cat b.cat ....\n\
\n\
DESCRIPTION\n\
	'pastecats' reads a set of catalogues and combines them\n\
	line by line into a single catalogue.\n\
	It reads the object contents information from each of\n\
	the input headers in turn, and installs each item into\n\
	the output catalogue's list of object items.  If no object\n\
	item names are repeated then the result is essentially the same\n\
	as using the unix 'paste' facility.  If an item name appears\n\
	more than once then the output catalogue will contain a single\n\
	object item of that name, whose value comes from the last input\n\
	object item of that name (though the position of the object\n\
	in the list of items is determined by the first occurrence).\n\
	The input catalogues will normally have the same length\n\
	and be in some ordered correspondence with each other.\n\
	Pastecats will return with exit value -1 if it appears\n\
	that the catalogues do not have identical numbers of objects.\n\
\n\
	pastecats uses the 'iostream' mechanism, so you can provide\n\
	catalog generating commands (terminated by the pipe symbol)\n\
	in place of catalog names.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@cita.utoronto.ca\n\
\n\n"


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include "../utils/error.h"
#include "../catlib/cat.h"
#include "../utils/iostream.h"


main(int argc, char *argv[])
{
	int		arg, ncats, i, nread, filetype;
	iostream	**ipstream;
	cathead		**ipcat, *opcat;
	object		**ipobj, *opobj;
	item		*ipitem, *opitem;

	ncats = argc - 1;
	if (!ncats || ((argv[1][0] == '-') && (strlen(argv[1]) > 1)))
		error_exit(usage);

	/* open the cats */
	ipstream = (iostream **) calloc(ncats, sizeof(iostream *));

	for (i = 0; i < ncats; i++) {
		if (!(ipstream[i] = openiostream(argv[i + 1], "r"))) {
			fprintf(stderr, "pastecats: Can't open %s\n", argv[i + 1]);
			exit(-1);
		}
	}

	/* read the cat headers */
	ipcat = (cathead **) calloc(ncats, sizeof(cathead *));
	ipobj = (object **) calloc(ncats, sizeof(object *));
	for (i = 0; i < ncats; i++){
		setcatipf((ipstream[i])->f);
		ipcat[i]  = readcathead();
	}

	/* create the output cathead */
	opcat = (cathead *) calloc(1, sizeof(cathead));
	copyheaderinfo(opcat, ipcat[0]);

	/* add comment */
	addargscomment(argc, argv, opcat);
	
	/* add the object items */
	for (i = 0; i < ncats; i++) {
		ipitem = ipcat[i]->objectitembase;
		while (ipitem) {
			allocitemcontents(ipitem, &(ipitem->addr), 0);
			opitem = copyitem(ipitem);
			opitem->addr = ipitem->addr;
			addobjectitem(opitem, opcat);
			ipitem = ipitem->next;
		}
		ipobj[i] = newobject(ipcat[i]);
		connectobjecttocathead(ipobj[i]);
	}

	/* and write it out */
	getcatipfiletype(&filetype);
	setcatopfiletype(filetype);
	writecathead(opcat);

	/* and create the opobj */
	opobj = newobject(opcat);
	connectobjecttocathead(opobj);

	while (1) {
		nread = 0;
		for (i = 0; i < ncats; i++) {
			setcatipf((ipstream[i])->f);
			if (readobject(ipobj[i]))
				nread++;
		}
		if (nread == ncats) {
			writeobject(opobj);
		} else {
			if (nread == 0) {
				for (i = 0; i < ncats; i++) {
					closeiostream(ipstream[i]);
				}
				exit(0);
			} else {
				exit(-1);
			}
		}
	}
}
