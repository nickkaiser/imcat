/*
 * lc.c main function for lc catalogue lister
 *
 * see lc.man for info
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../catlib/cat.h"
#include "getop.h"
#include "stack.h"
#include "operators.h"
#include "rpn.h"
#include "error.h"
#include "lc.h"

/* input modes */
#define STD_INPUT_MODE		0	/* read a standard format catalogue */
#define TAB_INPUT_MODE		1	/* convert a non-standard format table */

/* output modes */
#define STD_OUTPUT_MODE		0
#define COUNT_MODE	1
#define SORT_MODE	2
#define HEADER_ONLY_MODE	3
#define NO_HEADER_MODE		4

/* header modes */
#define FULL_HEADER_MODE	0
#define LABELS_ONLY_MODE	1

/* comment chopping mode */
#define	NO_CHOP_MODE 		0
#define CHOP_HEAD_MODE		1
#define	CHOP_TAIL_MODE		2

/* sort modes */
#define SORT_ASCENDING		0
#define SORT_DESCENDING		1

#define MAX_FUNCTIONS		1024

/*
#define printmanpage  system(MANPAGE); exit(0);
*/

static	cathead		*inputcathead, *outputcathead;
static	object		*inputobject = NULL, *outputobject;
static	int		nrpn = 0, copycomments, sortmode;
static	char		commenttag = '#';
static	double		inputobjectcount = 0, outputobjectcount = 0;
static	char		*inputcountname = NULL, *outputcountname = NULL;
static	rpnfunc		*rpnfunction[MAX_FUNCTIONS];
static	char		*rpnname[MAX_FUNCTIONS], *rpnstring[MAX_FUNCTIONS];

main(int argc, char *argv[]){
	int		arg, inputmode, outputmode;
	int		conditionaloutput, conditionarg, irpn, theindex, addhistory;
	double		*conditionval, thenumber, *sortval;
	rpnfunc		*condition, *sortfunction;
	item		*theitem;
	char		tempword[128], *theword, *sortstring, *imcatdir;
	int		commentchopmode, nchop, index;
	sortobject	*sortobjectbase = NULL, *thesortobject, **sortarray;
	void		**theaddrptr, *theaddr;
	int		headermode, inputfiletype, tempfiletype;
	long		seed;
	void		srand48(long seedval);

	/* defaults */
	inputmode = STD_INPUT_MODE;
	outputmode = STD_OUTPUT_MODE;
	conditionaloutput = 0;
	addhistory = 1;
	copycomments = 1;
	commentchopmode = NO_CHOP_MODE;
	headermode = FULL_HEADER_MODE;

	/* parse args */
	arg = 1;
	/* first we check to see if argv[1] is -u */
	if ((argc > 1) &&(!strncmp(argv[1], "-u", 2))) {
		imcatdir = getenv("IMCATDIR");
		sprintf(tempword, "cat %s/src/cattools/lc/lc.man 1>&2", imcatdir);
		system(tempword);
		exit(-1);
	}
	/* now we read the input cat header */
	/* first see if argv[1]  is -C option */
	if ((argc > 1) && (!strncmp(argv[1], "-C", 2))) {	/* converting unformatted table */
		inputmode = TAB_INPUT_MODE;
		inputcathead = (cathead *) calloc(1, sizeof(cathead));
		inputfiletype = TEXT_FILE_TYPE;
		arg++;
	} else {
		inputcathead = readcathead();
		getcatipfiletype(&inputfiletype);
	}

	/* now we create space for contents */
	for (index = 0; index < inputcathead->nobjectitems; index++) {
		allocitemcontents((inputcathead->itemlist)[index], &(((inputcathead->itemlist)[index])->addr), 0);
	}

 	/* this is for rpn function dereferencing value names */
	setsourcecathead(inputcathead);

	/* and create the output cat header and copy over header info from input cat */
	outputcathead = (cathead *) calloc(1, sizeof(cathead));
	copyheaderinfo(outputcathead, inputcathead);

	/* now process command line args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {	/* it must be a name or rpn assignment */
			parsename(argv[arg]);
			argv[arg] = quote(argv[arg]);
		} else {			/* its an option */
			switch (argv[arg][1]) {		/* big switch */
				case 'h':
					outputmode = HEADER_ONLY_MODE;
					break;
				case 'b':
					setcatopfiletype(BINARY_FILE_TYPE);
					break;
				case 'B':
					setcatopfiletype(inputfiletype);
					break;
				case 'x':
					addhistory = 0;
					break;
				case 'c':
					outputmode = COUNT_MODE;
					break;
				case 'o':
					outputmode = NO_HEADER_MODE;
					break;
				case 'O':
					headermode = LABELS_ONLY_MODE;
					break;
				case 'p':
					if (!(theitem = getheaderitem(argv[++arg], inputcathead))) {
						error_exit("lc: header item missing\n");
					}
					getcatopfiletype(&tempfiletype);
					setcatopfiletype(TEXT_FILE_TYPE);
					writeitem(theitem, theitem->addr, 0);
					setcatopfiletype(tempfiletype);
					fprintf(stdout, "\n");
					exit(0);
					break;
				case 'P':
					if (!(theitem = getheaderitem(argv[++arg], inputcathead))) {
						error_exit("lc: header item missing\n");
					}
					getcatopfiletype(&tempfiletype);
					setcatopfiletype(TEXT_FILE_TYPE);
					writeitem(theitem, theitem->addr, 0);
					setcatopfiletype(tempfiletype);
					fprintf(stdout, "\n");
					break;
				case 'q':
					sscanf(argv[++arg], "%ld", &seed);
					srand48(seed);
					break;
				case 'i':
					conditionaloutput = 1;
					arg++;
					conditionarg = arg;
					break;
				case 's':
					sortstring = argv[++arg];
					argv[arg] = quote(sortstring);
					outputmode = SORT_MODE;
					sortmode = SORT_ASCENDING;
					break;
				case 'S':
					sortstring = argv[++arg];
					argv[arg] = quote(sortstring);
					outputmode = SORT_MODE;
					sortmode = SORT_DESCENDING;
					break;
				case 'H':
					parseheaderexpr(argv[++arg]);
					argv[arg] = quote(argv[arg]);
					break;
				case 'r':
					deleteobjectitem(argv[++arg], outputcathead);
					break;
				case 'R':
					deleteheaderitem(argv[++arg], outputcathead);
					break;
				case 'a':
					addcomment(argv[++arg], outputcathead);
					break;
				case 'd':
					adddate();
					break;
				case 'L':
					commentchopmode = CHOP_HEAD_MODE;
					sscanf(argv[++arg], "%d", &nchop);
					nchop++;
					break;
				case 'F':
					commentchopmode = CHOP_TAIL_MODE;
					sscanf(argv[++arg], "%d", &nchop);
					break;
				case 'n':
					theitem = newitem(argv[++arg], NUM_TYPE, 1, 1);
					allocitemcontents(theitem, &(theitem->addr), 0);
					addobjectitem(theitem, inputcathead);
					theaddr = theitem->addr;
					theitem = copyitem(theitem);
					addobjectitem(theitem, outputcathead);
					theitem->addr = theaddr;
					break;
				case 't':
					theitem = newitem(argv[++arg], TEXT_TYPE, 1, 1);
					allocitemcontents(theitem, &(theitem->addr), 0);
					addobjectitem(theitem, inputcathead);
					theaddr = theitem->addr;
					theitem = copyitem(theitem);
					addobjectitem(theitem, outputcathead);
					theitem->addr = theaddr;
					break;
				case 'N':
					theitem = sscannewitem(NUM_TYPE, argv[++arg]);
					allocitemcontents(theitem, &(theitem->addr), 0);
					addobjectitem(theitem, inputcathead);
					theaddr = theitem->addr;
					theitem = copyitem(theitem);
					addobjectitem(theitem, outputcathead);
					theitem->addr = theaddr;
					argv[arg] = quote(argv[arg]);
					break;
				case 'T':
					theitem = sscannewitem(TEXT_TYPE, argv[++arg]);
					allocitemcontents(theitem, &(theitem->addr), 0);
					addobjectitem(theitem, inputcathead);
					theaddr = theitem->addr;
					theitem = copyitem(theitem);
					addobjectitem(theitem, outputcathead);
					theitem->addr = theaddr;
					argv[arg] = quote(argv[arg]);
					break;
				case 'Q':
					commenttag = argv[++arg][0];
					break;
				case 'I':
					copycomments = 0;
					break;
				case 'z':
					/* set swapbytesi for all the object items */
					theitem = inputcathead->objectitembase;
					while (theitem) {
						theitem->swapbytesi = 1;
						theitem = theitem->next;
					}
					break;
				default:
					fprintf(stderr, "lc: bad argument: do 'lc -u' to get man page\n");
					exit(-1);
					break;
			}
		}
		arg++;
	}
	
	if (inputmode == TAB_INPUT_MODE) {
		readtabheader();
	}

        inputobject = newobject(inputcathead);
        connectobjecttocathead(inputobject);

	if (conditionaloutput) {	
		condition = newrpnfunction("CONDITION", argv[conditionarg]);
		evalrpnfunction(condition);
		conditionval = (double *) ((condition->result)->addr);
		argv[conditionarg] = quote(argv[conditionarg]);
	}

	if (commentchopmode != NO_CHOP_MODE)
		chopcomments(outputcathead, commentchopmode, nchop);

	if (addhistory) {
		addargscomment(argc, argv, outputcathead);
	}

	/* if we have no objects in output cat we copy them over from input */
	if (!(outputcathead->nobjectitems)) {
		outputcathead->nobjectitems = inputcathead->nobjectitems;
		outputcathead->objectitembase = inputcathead->objectitembase;
		outputcathead->itemlist = inputcathead->itemlist;		
/*
		dospecialname("+all");
*/
	}



	switch (outputmode) {
		case STD_OUTPUT_MODE:
		case SORT_MODE:
			if (headermode == FULL_HEADER_MODE) {
				writecathead(outputcathead);
			} else {
				fprintf(stdout, "#");
				theitem = outputcathead->objectitembase;
        			while(theitem) {
					writelabel(theitem, 0);
					theitem = theitem->next;
				}
				fprintf(stdout, "\n");
			}
			break;
		case HEADER_ONLY_MODE:
			writecathead(inputcathead);
			exit(0);
			break;
		case COUNT_MODE:
		case NO_HEADER_MODE:
			break;
		default:
			error_exit("main: bad output mode\n");
			break;
	}
	
	/* create a new output object */
        outputobject = newobject(outputcathead);
	connectobjecttocathead(outputobject);

	if (outputmode == SORT_MODE) {
		/*setsourcecathead(outputcathead); */
		sortfunction = newrpnfunction("SORTFN", sortstring);
		evalrpnfunction(sortfunction);
		sortval = (double *) ((sortfunction->result)->addr);
	}

	while (readobject(inputobject)) {
		inputobjectcount++;
		if (conditionaloutput) {
			evalrpnfunction(condition);
			if (!((int) *conditionval))
				continue;
		}
		for (irpn = 0; irpn < nrpn; irpn++) {
			evalrpnfunction(rpnfunction[irpn]);
		}
		outputobjectcount++;
		switch (outputmode) {
			case STD_OUTPUT_MODE:
			case NO_HEADER_MODE:
				writeobject(outputobject);
				break;
			case SORT_MODE:
				thesortobject = (sortobject *) calloc(1, sizeof(sortobject));
				thesortobject->obj = newobject(outputcathead);
				for (index = 0; index < outputcathead->nobjectitems; index++) {
					theitem = (outputcathead->itemlist)[index];
					theaddrptr = (thesortobject->obj)->addrlist +index;
					allocitemcontents(theitem, theaddrptr, 0);
					copyitemcontents(theitem, *theaddrptr, theitem, (outputobject->addrlist)[index], 0);
				}
				evalrpnfunction(sortfunction);
				thesortobject->sortval = *sortval;			
				thesortobject->next = sortobjectbase;
				sortobjectbase = thesortobject;
				break;
			case COUNT_MODE:
				break;
		}	
	}

	
	if (outputmode == SORT_MODE) {
		/* create array */
		sortarray = (sortobject **) calloc(outputobjectcount, sizeof(sortobject *));
		/* load the sortobjects into the array */
		thesortobject = sortobjectbase;
		index = 0;
		while (thesortobject) {
			sortarray[index++] = thesortobject;
			thesortobject = thesortobject->next;
		}
		/* sort the array */
		qsort((char *) sortarray, outputobjectcount, sizeof(void *), objcmp);
		/* and output */
		for (index = 0; index < outputobjectcount; index++) {
			writeobject((sortarray[index])->obj);
		}
	}

	/* finally output count if in  COUNT_MODE */
	if (outputmode == COUNT_MODE)
		fprintf(stdout, "%d\n", (int) outputobjectcount);
	exit(0);
}




void	parsename(char *string)
{
	item	*inputitem, *outputitem;
	char	*stringbase, *word, *name;

	stringbase = string;
	getword(&string);
	if (!getword(&string)) {		/* simple named argument */
		if (stringbase[0] == '+') {
			dospecialname(stringbase);
		} else {
			inputitem = getobjectitem(stringbase, inputcathead);
			outputitem = copyitem(inputitem);
			addobjectitem(outputitem, outputcathead);
			outputitem->addr = inputitem->addr;
		}
	} else {				/* multi-word string */
		rpnname[nrpn] = getword(&stringbase);
		word = getword(&stringbase);
		if (strcmp(word, "=")) {
			error_exit("parsename: expected '='\n");
		}
		rpnstring[nrpn] = (char *) calloc(1 + strlen(stringbase), sizeof(char));
		strcpy(rpnstring[nrpn], stringbase);
		rpnfunction[nrpn] = newrpnfunction(rpnname[nrpn], rpnstring[nrpn]);
		evalrpnfunction(rpnfunction[nrpn]);
		addobjectitem((rpnfunction[nrpn])->result, outputcathead);
		if (++nrpn >= MAX_FUNCTIONS) {
			error_exit("parsename: too many functions\n");
		}
	}
}

void	parseheaderexpr(char *string)
{
	item		*theitem;
	char		*word, *name;
	rpnfunc 	*headerfunc;

	name = getword(&string);
	word = getword(&string);
	if (!word || (strcmp(word, "="))) {
		error_exit("parseheaderexpr: expected -H 'name = expr'");
	}
	headerfunc = newrpnfunction(name, string);
	evalrpnfunction(headerfunc);
	installitem(headerfunc->result, &(outputcathead->headeritembase));
	evalrpnfunction(headerfunc);
}	





void	addheadercontents(item *theitem, char *string, cathead *thecathead)
{
	allocitemcontents(theitem, &(theitem->addr), 0);
	sscanitem(theitem, theitem->addr, 0, &string);
	installitem(theitem, &(thecathead->headeritembase));
}

char	*quote(char *string)
{
	char	*result;
	
	result = (char *) calloc(strlen(string) + 3, sizeof(char));
	sprintf(result, "'%s'", string);
	return (result);
}



void	readtabheader(void)
{
	int	c;
	int	newline = 1, linepos = 0;
	char	line[1024];

	while (c = getc(stdin)) {
		if (newline && (c != commenttag)) {
			ungetc(c, stdin);
			return;
		} else {
			newline = 0;
		}
		if (c == '\n') {
			newline = 1;
			line[linepos++] = '\0';
			if (copycomments)
				addcomment(line, outputcathead);
			linepos = 0;
		} else {
			line[linepos++] = (char) c;
		}
	}
}



void	dospecialname(char *string)
{
	int	index;
	item	*theitem;

	string++;
	if (!strcmp(string, "all")) {
		for (index = 0; index < inputcathead->nobjectitems; index++) {
			theitem = copyitem((inputcathead->itemlist)[index]);
			addobjectitem(theitem, outputcathead);
			theitem->addr = ((inputcathead->itemlist)[index])->addr;
		}
	} else {
		theitem = newitem(string + 1, NUM_TYPE, 1, 1);
		addobjectitem(theitem, outputcathead);
		switch (string[0]) {
			case 'c':
				theitem->addr = &outputobjectcount;
				break;
			case 'C':
				theitem->addr = &inputobjectcount;
				break;
			default:
				error_exit("dospecialname: bad name\n");
				break;
		}
	}
}





void		chopcomments(cathead *thecathead, int mode, int nchop)
{
	int	count = 0;
	comment	*thecomment;

	if (nchop) {
		thecomment = thecathead->commentbase;
		while (thecomment) {
			count++;
			if (count >= nchop) {
				switch (mode) {
					case CHOP_TAIL_MODE:
						thecomment->next = NULL;
						break;
					case CHOP_HEAD_MODE:
						thecathead->commentbase = (thecathead->commentbase)->next;
						break;
					default:
						error_exit("chopcomments: bad chopping mode\n");
						break;
				}
			}
			thecomment = thecomment->next;
		}
	} else {
		thecathead->commentbase = NULL;
	}	
}




void	adddate(void){
	FILE		*datepipe;
	char		*datestring;	
	int	 	pos = 0;
	item		*theitem;	

	datestring = (char *) calloc(128, sizeof(char));
	if (!(datepipe = popen("date", "r"))) {
		error_exit("adddate: cannot open pipe to date command\n");
	}
	fgets(datestring, 124, datepipe);
	fclose(datepipe);
	while (1) {
		if (datestring[pos] == '\0')
			break;
		if (datestring[pos] == ' ')
			datestring[pos] = '_';		
		if (datestring[pos] == '\n')
			datestring[pos] = '\0';
		pos++;	
	}
	theitem = newitem("DATE", TEXT_TYPE, 1, 1);
	allocitemcontents(theitem, &(theitem->addr), 0);
	*((char **) (theitem->addr)) = datestring;
	installitem(theitem, &(outputcathead->headeritembase));
}


int	objcmp(const void *ptr1, const void *ptr2)
{
	sortobject 	**sortobj1, **sortobj2;

	sortobj1 = (sortobject **) ptr1;
	sortobj2 = (sortobject **) ptr2;
	if (sortmode == SORT_ASCENDING) {
		if ((*sortobj1)->sortval > (*sortobj2)->sortval) 
			return (1);
		if ((*sortobj1)->sortval < (*sortobj2)->sortval) 
			return (-1);
		return (0);
	} else {
		if ((*sortobj1)->sortval > (*sortobj2)->sortval) 
			return (-1);
		if ((*sortobj1)->sortval < (*sortobj2)->sortval) 
			return (1);
		return (0);
	}
}

