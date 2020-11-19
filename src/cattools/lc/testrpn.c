/*
 * rpn - test the stack functions
 */


#include <stdio.h>

#include "../../catlib/cat.h"
#include "getop.h"
#include "stack.h"
#include "operators.h"
#include "rpn.h"
#include "error.h"


main(int argc, char *argv[]){
        cathead *inputcathead, *outputcathead;
        object  *inputobject, *outputobject;
	item	*theitem;
	rpnfunc	*thefunction;
	int	theindex;
	char	*rpnstring, *word;
	
	if (argc != 2) {
		error_exit("usage: rpn string\n");
	}


        inputcathead = readcathead();
	setsourcecathead(inputcathead);
	inputobject = newobject(inputcathead);
        allocobjectcontents(inputobject);

	outputcathead = (cathead *) calloc(1, sizeof(cathead));
 	copyheaderinfo(outputcathead, inputcathead);
/*	deleteobjectitem("Q", outputcathead);*/

 	thefunction = newrpnfunction("RPNFUNC", argv[1]);
	evalrpnfunction(thefunction);
	addobjectitem(thefunction->result, outputcathead);

	theitem = newitem("FRED", TEXT_TYPE, 1, 1);
	addobjectitem(theitem, outputcathead);

	writecathead(outputcathead);

	outputobject = newobject(outputcathead);
        inheritcontents(outputobject, inputobject);
	theindex = getobjectitemindex("RPNFUNC", outputobject);
	setaddress(outputobject, theindex, (thefunction->result)->addr);
        allocobjectcontents(outputobject);
        while (readobject(inputobject)) {
 		evalrpnfunction(thefunction);
                writeobject(outputobject);	
        }

}
