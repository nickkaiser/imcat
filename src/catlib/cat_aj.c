/*
 * cat.c
 *
 * catalogue I/O functions
 *
 */



#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "cat.h"
#include "error.h"

/* 
   AOJ  3 jun 1999 -- moved initialization of catipf and catopf
                      from here to each function using the pnts
*/
static FILE	*catipf;
static FILE	*catopf;

static int	inputfiletype = TEXT_FILE_TYPE;
static int	outputfiletype = TEXT_FILE_TYPE;

static int	swapbytesi = 0;
static int	swapbyteso = 0;


#define N_TYPES		8
static	char	*type[N_TYPES] = {"number", "text", "comment:", "contents:",
"op2", "op1", "genericop", "rpn_function"};


#define	WORD_SIZE		128
#define	LINE_SIZE		100000

#define LINE1_TEXT	"# text format catalogue file --- do not edit manually, use 'lc'\n"
#define LINE1_BINARY	"# binary catalogue file --- use 'lc' to list contents\n"
#define LINE2	"# header:\n"


cathead		*readcathead(void)
{
  item	        *theitem;
  cathead	*thecathead;
  int	        itype, done = 0, i, filetype;
  char	        line[LINE_SIZE];
  comment	*thecomment;

  /* AOJ 990615 -- initialize catipf if not already done */
  if (!catipf)
    catipf = stdin;

  thecathead = (cathead *) calloc(1, sizeof(cathead));
  /* check the first 2 lines */
  fgets(line, LINE_SIZE, catipf);
  filetype = TEXT_FILE_TYPE;
  setcatipfiletype(TEXT_FILE_TYPE);
  if (strcmp(line, LINE1_TEXT)) {
    filetype = BINARY_FILE_TYPE;
    if (strcmp(line, LINE1_BINARY)) {
      error_exit("readcathead: file format error\n");
		}
  }
  fgets(line, LINE_SIZE, catipf);
  if (strcmp(line, LINE2)) {
		error_exit("readcathead: file format error\n");
  }
  while (!done) {
    switch (itype = gettype()) {
    case NUM_TYPE:
    case TEXT_TYPE:
      theitem = readnewitem(itype);
      allocitemcontents(theitem, &(theitem->addr), 0);
      readitem(theitem, theitem->addr, 0);
      installitem(theitem, &(thecathead->headeritembase));
      break;
    case COMMENT_TYPE:
      fgets(line, LINE_SIZE, catipf);
      thecomment = (comment *) calloc(1, sizeof(comment));
      thecomment->text = (char *) calloc(strlen(line)+1, sizeof(char));
      strcpy(thecomment->text, line);
      installcomment(thecomment, &(thecathead->commentbase));
      break;
    case CONTENTS_TYPE:
      fscanf(catipf, "%d", &(thecathead->nobjectitems));
      for (i = 0; i < thecathead->nobjectitems; i++) {
	itype = gettype();
	theitem = readnewitem(itype);
	installitem(theitem, &(thecathead->objectitembase));
      }
      makeitemlist(thecathead);
      done = 1;
      break;
    default:
      fprintf(stderr, "readnewitem: bad type\n");
      exit(-1);
      break;
    }
  }
  fgets(line, LINE_SIZE, catipf);
  fgets(line, LINE_SIZE, catipf);
  setcatipfiletype(filetype);
  return(thecathead);
}


void		copyheaderinfo(cathead *thecathead, cathead *sourcecathead)
{
  item	*theitem, *sourceitem;
  comment	*thecomment, *sourcecomment;

  sourceitem = sourcecathead->headeritembase;
  while (sourceitem) {
    theitem = copyitem(sourceitem);
    allocitemcontents(theitem, &(theitem->addr), 0);
    copyitemcontents(theitem, theitem->addr, sourceitem, sourceitem->addr, 0);
    installitem(theitem, &(thecathead->headeritembase));
    sourceitem = sourceitem->next;
  }
  sourcecomment = sourcecathead->commentbase;
  while (sourcecomment) {
    thecomment = copycomment(sourcecomment);
    installcomment(thecomment, &(thecathead->commentbase));
    sourcecomment = sourcecomment->next;
  }
}



void		copycontentinfo(cathead *thecathead, cathead *sourcecathead)
{
  item	*theitem, *sourceitem;
  int	index = 0;
  
  thecathead->nobjectitems = sourcecathead->nobjectitems;
  thecathead->itemlist = (item **) calloc(thecathead->nobjectitems, sizeof(item *));
  sourceitem = sourcecathead->objectitembase;
  while (sourceitem) {
    theitem = copyitem(sourceitem);
    allocitemcontents(theitem, &(theitem->addr), 0);
    copyitemcontents(theitem, theitem->addr, sourceitem, sourceitem->addr, 0);
    installitem(theitem, &(thecathead->objectitembase));
    (thecathead->itemlist)[index++] = theitem;
    sourceitem = sourceitem->next;
  }	
}

void		writecathead(cathead *thecathead)
{
  item	*theitem;
  comment	*thecomment;
  int	filetype;

  /* AOJ 990615 -- initialize catopf if not already done */
  if (!catopf)
    catopf = stdout;
  
  getcatopfiletype(&filetype);
  setcatopfiletype(TEXT_FILE_TYPE);
  switch (filetype) {
  case TEXT_FILE_TYPE:
    fprintf(catopf, LINE1_TEXT);
    break;
  case BINARY_FILE_TYPE:
    fprintf(catopf, LINE1_BINARY);
    break;
  default:
    error_exit("writecathead: bad output file type\n");
    break;
  }
  fprintf(catopf, LINE2);
  theitem = thecathead->headeritembase;
  while(theitem) {
    writename(theitem);
    writeitem(theitem, theitem->addr, 0);
    fprintf(catopf, "\n");
    theitem = theitem->next;
  }
  thecomment = thecathead->commentbase;
  while (thecomment) {
    fprintf(catopf, "# comment:%s", thecomment->text);
    thecomment = thecomment->next;
  }
  fprintf(catopf, "# contents: %d\n", thecathead->nobjectitems);
  theitem = thecathead->objectitembase;
  while(theitem) {
    writename(theitem);
    fprintf(catopf, "\n");
    theitem = theitem->next;
  }
  fprintf(catopf, "#");
  theitem = thecathead->objectitembase;
  while(theitem) {
    writelabel(theitem, 0);
    theitem = theitem->next;
  }
  fprintf(catopf, "\n");
  setcatopfiletype(filetype);
}


void		freecathead(cathead *thecathead)
{
  item	*theitem, *nextitem;
  comment	*thecomment, *nextcomment;

  /* trash the header items */
  theitem = thecathead->headeritembase;
  while (theitem) {
    nextitem = theitem->next;
    freeitem(theitem);
    theitem = nextitem;
  }
  /* trash the object items */
  theitem = thecathead->objectitembase;
  while (theitem) {
    nextitem = theitem->next;
    freeitem(theitem);
    theitem = nextitem;
  }
  free(thecathead->itemlist);
}

void		makeitemlist(cathead *thecathead)
{
  item	*theitem;
  int	index = 0;
  
  if (thecathead->itemlist)
    free(thecathead->itemlist);
  thecathead->itemlist = (item **) calloc(thecathead->nobjectitems, sizeof(item *));
  theitem = thecathead->objectitembase;
  while (theitem) {
    (thecathead->itemlist)[index] = theitem;
    theitem = theitem->next;
    index++;
  } 
}


int		gettype(void)
{
  int	itype;
  char	word[WORD_SIZE];

  /* AOJ 990615 -- initialize catipf if not already done */
  if (!catipf)
    catipf = stdin;

  /* ignore # character and fetch next word */
  while (1) {
    fscanf(catipf, "%s", word);
    /* printf("debug: word = %s\n", word);*/
    if (strcmp(word, "#"))
      break;
  }
  for (itype = 0; itype < N_TYPES; itype++) {
    if (!strcmp(word, type[itype]))
      break;
  }
  if (itype == N_TYPES) {
    fprintf(stderr, "gettype: input format error\n");
    exit(-1);
  }
  return (itype);
}


item		*readnewitem(int itype)
{
  item	*theitem;
  int	ndim, idim, dim;
  char	name[WORD_SIZE];

  /* AOJ 990615 -- initialize catipf if not already done */
  if (!catipf)
    catipf = stdin;
  
  theitem = (item *) calloc(1, sizeof(item));
  theitem->itype = itype;
  fscanf(catipf, "%d", &ndim);
  theitem->ndim = ndim;
  for (idim = 0; idim < ndim; idim++) {
    if (idim >= MAX_DIMS) {
      fprintf(stderr, "readnewitem: too many dimensions\n");
      exit(-1);
    }
    fscanf(catipf, "%d", &dim);
    theitem->dim[idim] = dim;
  }
  fscanf(catipf, "%s", name);
  theitem->name = (char *) calloc(strlen(name) + 1, sizeof(char));
  strcpy(theitem->name, name);
  return (theitem);
}

item		*sscannewitem(int itype, char *string)
{
  item	*theitem;
  int	ndim, idim, dim, shift;
  char	name[WORD_SIZE];

  theitem = (item *) calloc(1, sizeof(item));
  theitem->itype = itype;
  sscanf(string, "%d%n", &ndim, &shift);
  string += shift;
  theitem->ndim = ndim;
  for (idim = 0; idim < ndim; idim++) {
    if (idim >= MAX_DIMS) {
      fprintf(stderr, "readnewitem: too many dimensions\n");
      exit(-1);
    }
    sscanf(string, "%d%n", &dim, &shift);
    string += shift;
    theitem->dim[idim] = dim;
  }
  sscanf(string, "%s%n", name, &shift);
  string += shift;
  theitem->name = (char *) calloc(strlen(name) + 1, sizeof(char));
  strcpy(theitem->name, name);
  return (theitem);
}


item		*newitem(char *name, int itype, int ndim, ...)
{
  va_list	ap;
  item	*theitem;
  int	idim;
  
  if ((itype < 0) || (itype >= N_TYPES)) {
    fprintf(stderr, "newitem: illegal type requested\n");
    exit(-1);
  } 
  theitem = (item *) calloc(1, sizeof(item));
  theitem->name = (char *) calloc(1 + strlen(name), sizeof(char *));
  strcpy(theitem->name, name);
  theitem->itype = itype;
  theitem->ndim = ndim;
  va_start(ap, ndim);
  for (idim = 0; idim < ndim; idim++) {
    (theitem->dim)[idim] = va_arg(ap, int);
  }
  return (theitem);
}


item		*newitembydimarray(char *name, int itype, int ndim, int *dim)
{
  item	*theitem;
  int	idim;
  
  if ((itype < 0) || (itype >= N_TYPES)) {
    fprintf(stderr, "newitem: illegal type requested\n");
    exit(-1);
  } 
  theitem = (item *) calloc(1, sizeof(item));
  theitem->name = (char *) calloc(1 + strlen(name), sizeof(char *));
  strcpy(theitem->name, name);
  theitem->itype = itype;
  theitem->ndim = ndim;
  for (idim = 0; idim < ndim; idim++) {
    (theitem->dim)[idim] = dim[idim];
  }
  return (theitem);
}



item		*copyitem(item *sourceitem)
{
  item	*theitem;
  
  theitem = (item *) calloc(1, sizeof(item));
  *theitem = *sourceitem;
  theitem->next = NULL;
  theitem->addr = NULL;
  if (sourceitem->name) {
    theitem->name = (char *) calloc(strlen(sourceitem->name) + 1, sizeof(char));
    strcpy(theitem->name, sourceitem->name);
  }
  return (theitem);
}



int		allocitemcontents(item *theitem, void **addr, int level)
{
  int	i;
  void	*theaddr;
  
  if (level < (theitem->ndim - 1)) {
    theaddr = (*addr) = (void *) calloc((theitem->dim)[level], sizeof(void *));
    for (i = 0; i < (theitem->dim)[level]; i++) {
      allocitemcontents(theitem, (void **) theaddr + i, level + 1);
    }
  } else {
    theaddr = *addr = (void *) calloc((theitem->dim)[level], size(theitem->itype));
  }
}


void		copyitemcontents(item *dstitem, void *dstaddr, item *srcitem, void *srcaddr, int level)
{
  int	i;
  char	*srcword;
  
  if (level < (srcitem->ndim - 1)) {
    for (i = 0; i < (srcitem->dim)[level]; i++) {
      copyitemcontents(dstitem, *((void **) dstaddr + i), srcitem, *((void **) srcaddr + i), level + 1);
    }
  } else {
    for (i = 0; i < (srcitem->dim)[level]; i++) {
      switch(dstitem->itype) {
      case NUM_TYPE:
	*((double *) dstaddr + i) = *((double *) srcaddr + i);
	break;
      case TEXT_TYPE:
	srcword = *((char **) srcaddr + i);
	if (srcword) {
	  *((char **) dstaddr + i) = (char *) calloc(1 + strlen(srcword), sizeof(char));
	  strcpy(*((char **) dstaddr + i), srcword);
	}
	break;
      default:
	fprintf(stderr, "copyitemcontents: bad type\n");
	exit(-1);
	break;
      }
    }
  }
}

void		freeitemcontents(item *theitem, void *addr, int level)
{
  int	i;
  
  if (level < (theitem->ndim - 1)) {
    for (i = 0; i < (theitem->dim)[level]; i++) {
      freeitemcontents(theitem, *((void **) addr + i), level + 1);
      free(addr);			
    }
  } else {
    free(addr);
  }
}

void		freeitem(item *theitem)
{
  if (theitem->name)
    free(theitem->name);
  theitem->name = NULL;
  if ( theitem->addr)
    freeitemcontents(theitem, theitem->addr, 0);
  free(theitem);
}


int		readitem(item *theitem, void *addr, int level)
{
  int	i, nmatch, len;
  char	word[WORD_SIZE], *slot;

  /* AOJ 990615 -- initialize catipf if not already done */
  if (!catipf)
    catipf = stdin;
  
  if (level < (theitem->ndim - 1)) {
    for (i = 0; i < (theitem->dim)[level]; i++) {
      if (!readitem(theitem, *((void **) addr + i), level + 1))
	return(0);
    }
  } else {
    switch (inputfiletype) {
    case TEXT_FILE_TYPE:
      for (i = 0; i < (theitem->dim)[level]; i++) {
	nmatch = fscanf(catipf, "%s", word);
	if (nmatch == 0 || nmatch == EOF)
	  return (0);
	switch (theitem->itype) {
	case NUM_TYPE:
	  if (1 != sscanf(word, "%lf", ((double *) addr + i))) {
	    error_exit("readitem: unable to convert input\n");
	  }
	  break;
	case TEXT_TYPE:
	  if (*((char **) addr + i))
	    free(*((char **) addr + i));
	  slot = *((char **) addr + i) = (char *) calloc(strlen(word) + 1, sizeof(char));
	  strcpy(slot, word);
	  break;
	default:
	  fprintf(stderr, "readitem: bad type\n");
	  exit(-1);
	  break;
	}
      }
      break;
    case BINARY_FILE_TYPE:
      switch (theitem->itype) {
      case NUM_TYPE:
	if ((theitem->dim)[level] != fread(addr, sizeof(double), (theitem->dim)[level], catipf))
	  return(0);
	if (swapbytesi) {
	  swapbytes(addr, (theitem->dim)[level]);
	}
	break;
      case TEXT_TYPE:
	for (i = 0; i < (theitem->dim)[level]; i++) {
	  slot = *((char **) addr + i);
	  if (slot)
	    free(slot);
	  if (1 != fread(&len, sizeof(int), 1, catipf))
	    return(0);
	  slot = *((char **) addr + i) = (char *) calloc(len + 1, sizeof(char));
	  if (1 != fread(slot, len * sizeof(char), 1, catipf))
	    return(0);
	}
	break;
      default:
	fprintf(stderr, "readitem: bad type\n");
	exit(-1);
	break;
      }
      break;
    default:
      error_exit("readitem: bad input file type\n");
      break;
    }
  }
  return (1);
}





int		sscanitem(item *theitem, void *addr, int level, char **string)
{
  int	i, nmatch, shift;
  char	word[WORD_SIZE], *slot;

  if (level < (theitem->ndim - 1)) {
    for (i = 0; i < (theitem->dim)[level]; i++) {
      if (!sscanitem(theitem, *((void **) addr + i), level + 1, string))
	return(0);
    }
  } else {
    for (i = 0; i < (theitem->dim)[level]; i++) {
      nmatch = sscanf(*string, "%s%n", word, &shift);
      *string += shift;
      if (nmatch == 0 || nmatch == EOF)
	return (0);
      switch (theitem->itype) {
      case NUM_TYPE:
	sscanf(word, "%lf", ((double *) addr + i));
	break;
      case TEXT_TYPE:
	if (*((char **) addr + i))
	  free(*((char **) addr + i));
	slot = *((char **) addr + i) = (char *) calloc(strlen(word) + 1, sizeof(char));
	strcpy(slot, word);
	break;
      default:
	fprintf(stderr, "readitem: bad type\n");
	exit(-1);
	break;
      }
    }
  }
  return (1);
}





int		writeitem(item *theitem, void *addr, int level)
{
  int	i, len;
  char	*slot;

  /* AOJ 990615 -- initialize catopf if not already done */
  if (!catopf)
    catopf = stdout;
  
  if (level < (theitem->ndim - 1)) {
    for (i = 0; i < (theitem->dim)[level]; i++) {
      writeitem(theitem, *((void **) addr + i), level + 1);
    }
  } else {
    switch (outputfiletype) {
    case TEXT_FILE_TYPE:
      for (i = 0; i < (theitem->dim)[level]; i++) {
	switch (theitem->itype) {
	case NUM_TYPE:
	  fprintf(catopf, LC_NUM_FMT, *((double *) addr + i));
	  break;
	case TEXT_TYPE:
	  slot = *((char **) addr + i);
	  if (!slot) {
	    fprintf(catopf, LC_TEXT_FMT, "(null)");
	  } else {
	    fprintf(catopf,  LC_TEXT_FMT, slot);
	  }
	  break;
	default:
	  fprintf(stderr, "writeitem: bad type\n");
	  exit(-1);
	  break;
	}
      }
      break;
    case BINARY_FILE_TYPE:
      switch (theitem->itype) {
      case NUM_TYPE:
	if (swapbyteso) {
	  swapbytes(addr, (theitem->dim)[level]);
	}
	fwrite(addr, sizeof(double), (theitem->dim)[level], catopf);
	break;
      case TEXT_TYPE:
	for (i = 0; i < (theitem->dim)[level]; i++) {
	  slot = *((char **) addr + i);
	  len = strlen(slot);
	  fwrite(&len, sizeof(int), 1, catopf);
	  fwrite(slot, len * sizeof(char), 1, catopf);
	}
	break;
      default:
	fprintf(stderr, "writeitem: bad type\n");
	exit(-1);
	break;
      }
      break;
    default:
      error_exit("writeitem: bad output file type\n");
      break;
    }
  }
}


int		writelabel(item *theitem, int level)
{
  int	i, dim;
  char	label[WORD_SIZE], word[WORD_SIZE];

  /* AOJ 990615 -- initialize catopf if not already done */
  if (!catopf)
    catopf = stdout;
  
  for (i = 0; i < (theitem->dim)[level]; i++) {
    (theitem->idim)[level] = i;
    if (level < (theitem->ndim - 1)) {
      writelabel(theitem, level + 1);
    } else {
      sprintf(label, "%s", theitem->name);
      if ((theitem->dim)[0] > 1) {
	for (dim = 0; dim < theitem->ndim; dim++) {
	  sprintf(word, "[%d]", (theitem->idim)[dim]);
	  strcat(label, word);
	}
      }
      if (theitem->itype == TEXT_TYPE) {
	fprintf(catopf, LC_TEXT_FMT, label);
      } else {
	fprintf(catopf, LC_RIGHT_TEXT_FMT, label);
      }	
    }
  }
}

int		writename(item *theitem)
{
  char	*theword, *thetype;
  int	i;

  /* AOJ 990615 -- initialize catopf if not already done */
  if (!catopf)
    catopf = stdout;
  
  theword = (char *) calloc(WORD_SIZE, sizeof(char));
  thetype = (char *) calloc(WORD_SIZE, sizeof(char));
  sprintf(thetype, "%-8s", type[theitem->itype]);
  sprintf(theword, " %d", theitem->ndim);
  strcat(thetype, theword);
  for (i = 0; i < theitem->ndim; i++) {
    sprintf(theword, " %d", (theitem->dim)[i]);
    strcat(thetype, theword);
  }
  fprintf(catopf, LC_HEADER_NAME_FMT, thetype, theitem->name);	
  free(theword);
  free(thetype);	
}


int		size(int itype)
{
  switch (itype) {
  case NUM_TYPE:
    return(sizeof(double));
    break;
  case TEXT_TYPE:
    return(sizeof(char *));
    break;
  default:
    fprintf(stderr, "size: bad type\n");
    exit(-1);
    break;
  }
}


int		installitem(item *theitem, item **base)
{
  item	*anitem, *last = NULL;
  int	replaced = 0;

  anitem = *base;
  if (!anitem) {
    *base = theitem;
  } else {
    while(anitem) {
      if (!strcmp(anitem->name, theitem->name)) {
	if (last) {
	  last->next = theitem;
	} else {
	  *base = theitem;
	}
	theitem->next = anitem->next;
	replaced = 1;
	/* freeitem(anitem);*/
	break;
      }
      last = anitem;
      anitem = anitem->next;
    }
    if (!replaced) {
      last->next = theitem;
    }
  }
  return (replaced);
}


void		installcomment(comment *thecomment, comment **base)
{
  comment	*acomment;

  acomment = *base;
  if (!acomment) {
    *base = thecomment;
  } else {
    while(acomment->next) {
      acomment = acomment->next;
    }
    acomment->next = thecomment;
  }
}


void		addcomment(char *thetext, cathead *thecathead)
{
  comment	*thecomment;

  thecomment = (comment *) calloc(1, sizeof(comment));
  thecomment->text = (char *) calloc(strlen(thetext) + 3, sizeof(char));
  strcpy(thecomment->text, " ");
  strcat(thecomment->text, thetext);
  strcat(thecomment->text, "\n");
  installcomment(thecomment, &(thecathead->commentbase));
}



void		addargscomment(int argc, char *argv[], cathead *thecathead)
{
  comment	*thecomment;
  int	len = 2, arg;

  thecomment = (comment *) calloc(1, sizeof(comment));
  for (arg = 0; arg < argc; arg++) {
    len += 1 + strlen(argv[arg]);
  }
  len += strlen(" history:");
  thecomment->text = (char *) calloc(len, sizeof(char));
  strcat(thecomment->text, " history:");
  for (arg = 0; arg < argc; arg++) {
    strcat(thecomment->text, " ");
    strcat(thecomment->text, argv[arg]);
  }
  strcat(thecomment->text, "\n");
  installcomment(thecomment, &(thecathead->commentbase));
}

comment		*copycomment(comment *sourcecomment)
{
  comment	*thecomment;

  thecomment = (comment *) calloc(1, sizeof(comment));
  if (sourcecomment->text) {
    thecomment->text = (char *) calloc(strlen(sourcecomment->text) + 1, sizeof(char));
    strcpy(thecomment->text, sourcecomment->text);
  }
  return (thecomment);
}


object		*newobject(cathead *thecathead)
{
  int	index;
  object	*theobject;

  theobject = (object *) calloc(1, sizeof(object));
  theobject->cathead = thecathead;
  theobject->nitems = thecathead->nobjectitems;
  theobject->addrlist = (void **) calloc(theobject->nitems, sizeof(void *));
  /*
    for (index = 0; index <= theobject->nitems; index++) {
    (theobject->addrlist)[index] = ((thecathead->itemlist)[index])->addr;
    }
  */
  return(theobject);
}


void		allocobjectcontents(object *theobject)
{
  int	index;
  item	*theitem;

  for (index = 0; index < theobject->nitems; index++) {
    if (!((theobject->addrlist)[index])) {
      theitem = ((theobject->cathead)->itemlist)[index];
      allocitemcontents(theitem, &((theobject->addrlist)[index]), 0);
      theitem->addr = (theobject->addrlist)[index];
    }
  }
}


void		connectcatheadtoobject(object *theobject)
{
  int	index;
  item	*theitem;

  for (index = 0; index < theobject->nitems; index++) {
    theitem = ((theobject->cathead)->itemlist)[index];
    theitem->addr = (theobject->addrlist)[index];
  }
}


void		connectobjecttocathead(object *theobject)
{
  int	index;
  item	*theitem;

  for (index = 0; index < theobject->nitems; index++) {
    theitem = ((theobject->cathead)->itemlist)[index];
    (theobject->addrlist)[index] = theitem->addr;
  }
}



int		readobject(object *theobject)
{
  int	index;
  item	*theitem;

  for (index = 0; index < theobject->nitems; index++) {
    theitem = ((theobject->cathead)->itemlist)[index];
    if (!readitem(theitem, (theobject->addrlist)[index], 0))
      return(0);
  }
  return(1);
}


void		writeobject(object *theobject)
{
  int	index;
  item	*theitem;

  /* AOJ 990615 -- initialize catopf if not already done */
  if (!catopf)
    catopf = stdout;
  
  if (outputfiletype == TEXT_FILE_TYPE)
    fprintf(catopf, " ");
  for (index = 0; index < theobject->nitems; index++) {
    theitem = ((theobject->cathead)->itemlist)[index];
    writeitem(theitem, (theobject->addrlist)[index], 0);
  }
  if (outputfiletype == TEXT_FILE_TYPE)
    fprintf(catopf, "\n");
}


void		freeobject(object *theobject)
{
  int		index;
  item	*theitem;
  
  for (index = 0; index < theobject->nitems; index++) {
    theitem = ((theobject->cathead)->itemlist)[index];
    freeitemcontents(theitem, &((theobject->addrlist)[index]), 0);
  }
  free(theobject->addrlist);
  free(theobject);
}



item		*getheaderitem(char *name, cathead *thecathead)
{
  item	*theitem;

  theitem = thecathead->headeritembase;
  while (theitem) {
    if (!strcmp(theitem->name, name)) {
      break;
    } else {
      theitem = theitem->next;
    }
    /*
    if (!theitem) {
      fprintf(stderr, "getheaderitem: item \"%s\" nonexistent\n", name);
      exit(-1);
    }
    */
  }
  return (theitem);
}



void		*getheaderitemaddress(char *name, cathead *thecathead)
{
	item	*theitem;

	theitem = getheaderitem(name, thecathead);
	if (theitem) {
		return(theitem->addr);
	} else {
		fprintf(stderr, "getheaderitem: item \"%s\" nonexistent\n", name);
		exit(-1);
	}
}


void		setheaderitemaddress(char *name, cathead *thecathead, void *theaddr)
{
  getheaderitem(name, thecathead)->addr = theaddr;
}




void		deleteheaderitem(char *name, cathead *thecathead)
{
  item	*theitem, *lastitem;
  
  lastitem = NULL;
  theitem = thecathead->headeritembase;
  while (theitem) {
    if (!strcmp(theitem->name, name)) {
      break;
    } else {
      lastitem = theitem;
      theitem = theitem->next;
    }
    if (!theitem) {
      fprintf(stderr, "getheaderitem: item \"%s\" nonexistent\n", name);
      exit(-1);
    }
  }
  if (lastitem) {
    lastitem->next = theitem->next;
  } else {
    thecathead->headeritembase = theitem->next;
  }
  freeitem(theitem);
}




void		addobjectitem(item *newitem, cathead *thecathead)
{
  if (!installitem(newitem, &(thecathead->objectitembase)))
    (thecathead->nobjectitems)++;
  makeitemlist(thecathead);		
}



void		deleteobjectitem(char *name, cathead *thecathead)
{
  item	*theitem, *lastitem;

  lastitem = NULL;
  theitem = thecathead->objectitembase;
  while (theitem) {
    if (!strcmp(theitem->name, name)) {
      break;
    } else {
      lastitem = theitem;
      theitem = theitem->next;
    }
    if (!theitem) {
      fprintf(stderr, "getheaderitem: item \"%s\" nonexistent\n", name);
      exit(-1);
    }
  }
  if (lastitem) {
    lastitem->next = theitem->next;
  } else {
    thecathead->objectitembase = theitem->next;
  }
  /* 	freeitem(theitem);*/
  thecathead->nobjectitems--;
  makeitemlist(thecathead);
}


item		*getobjectitem(char *name, cathead *thecathead)
{
  item	*theitem;

  theitem = thecathead->objectitembase;
  while (theitem) {
    if (!strcmp(theitem->name, name)) {
      break;
    } else {
      theitem = theitem->next;
    }
    if (!theitem) {
      fprintf(stderr, "getobjectitem: item \"%s\" nonexistent\n", name);
      exit(-1);
    }
  }
  return (theitem);
}



int		getobjectitemindex(char *name, object *theobject)
{
  item	*theitem;
  int	index = 0;

  for (index = 0; index < theobject->nitems; index++) {
    theitem = ((theobject->cathead)->itemlist)[index];
    if (!strcmp(theitem->name, name)) {
      break;
    }
  }
  if (index == theobject->nitems) {
    fprintf(stderr, "getobjectitemindex: item \"%s\" nonexistent\n", name);
    exit(-1);
  }
  return (index);
}




int		inheritcontents(object *theobject, object *sourceobject)
{
  int 	theindex, sourceindex, idim, identical = 1;
  char	*thename, *sourcename;
  item	*theitem, *sourceitem;
  
  for (theindex = 0; theindex < theobject->nitems; theindex++) {
    if (!((theobject->addrlist)[theindex])) {
      thename = (((theobject->cathead)->itemlist)[theindex])->name;
      for (sourceindex = 0; sourceindex < sourceobject->nitems; sourceindex++) {
	sourcename = (((sourceobject->cathead)->itemlist)[sourceindex])->name;
	if (!strcmp(thename, sourcename)) {
	  theitem = ((theobject->cathead)->itemlist)[theindex];
	  sourceitem = ((sourceobject->cathead)->itemlist)[sourceindex];
	  if (    theitem->itype != sourceitem->itype ||
		  theitem->ndim != sourceitem->ndim) {
	    identical = 0;
	  }
	  for (idim = 0; idim < theitem->ndim; idim++) {
	    if ((theitem->dim)[idim] != (sourceitem->dim)[idim]) {
	      identical = 0;
	    }
	  }
	  if (!identical) {
	    fprintf(stderr, "inheritcontents: non-identical contents type\n");
	    exit(-1);
	  }
	  (theobject->addrlist)[theindex] = (sourceobject->addrlist)[sourceindex];
	  break;
	}
      }
    }
  }
}




void		*getaddress(object *theobject, int theindex)
{
  if (theindex < 0 || theindex >= theobject->nitems) {
    fprintf(stderr, "getaddress: index out of allowed range\n");
    exit(-1);
  }
  return ((theobject->addrlist)[theindex]);
}




void		setaddress(object *theobject, int theindex, void *theaddr)
{
  if (theindex < 0 || theindex >= theobject->nitems) {
    fprintf(stderr, "getaddress: index out of allowed range\n");
    exit(-1);
  }
  (theobject->addrlist)[theindex] = theaddr;
}




void		copystring(char **thestring, char *sourcestring)
{
  if (*thestring)
    free(*thestring);
  *thestring = (char *) calloc(1 + strlen(sourcestring), sizeof(char));
  strcpy(*thestring, sourcestring);
}


void		setcatipf(FILE *stream)
{
  catipf = stream;
}


void		setcatopf(FILE *stream)
{
  catopf = stream;
}



void		getcatipfiletype(int *filetype)
{
  *filetype = inputfiletype;
}


void		setcatipfiletype(int filetype)
{
  inputfiletype = filetype;
}



void		getcatopfiletype(int *filetype)
{
  *filetype = outputfiletype;
}



void		setcatopfiletype(int filetype)
{
  outputfiletype = filetype;
}



void		swapbytes(void *addr, int ndoubles)
{
  static char	tmpchar[8];
  int		i, b;
  char		*caddr;

  caddr = (char *) addr;
  for (i = 0; i < ndoubles; i++) {
    for (b = 0 ; b < 8; b++) {
      tmpchar[b] = caddr[b];
    }
    for (b = 0 ; b < 8; b++) {
      caddr[b] = tmpchar[7 - b];
    }
    caddr += 8;
  }
}

void		swapbytesonoutput(void)
{
  swapbyteso = 1;
}


void		swapbytesoninput(void)
{
  swapbytesi = 1;
}


int	install1numericheaderitem(char *name, double value, cathead *thecat)
{
	item 	*theitem;

	theitem = newitem(name, NUM_TYPE, 1, 1);
	allocitemcontents(theitem, &(theitem->addr), 0);
	((double *)(theitem->addr))[0] = value;
	installitem(theitem, &(thecat->headeritembase));
}


int	install1textheaderitem(char *name, char *value, cathead *thecat)
{
	item 	*theitem;

	theitem = newitem(name, TEXT_TYPE, 1, 1);
	allocitemcontents(theitem, &(theitem->addr), 0);
	((char **)(theitem->addr))[0] = value;
	installitem(theitem, &(thecat->headeritembase));
}

