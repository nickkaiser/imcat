/*
 * fits.c
 *
 * functions to implement fits I/O
 */

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include "fits.h"
#include "convertarray.h"
#include "../utils/error.h"
#include "../utils/arrays.h"
#include "nan.h"


/* readfitsheader() reads a fits header from stream.
 *
 * intpixtype is set to FLOAT_PIXTYPE by default
 *
 */
fitsheader *readfitsheader(FILE *stream)
{
	fitsheader 	*fitshead;
	fitscomment	*newcomment, *thecomment;
	char		record[FITS_REC_SIZE], string[81], *textval, *byteswapenvvar;
	int		line, recorddone, alldone = 0, bad, dim, oldimcatformat;

	/* allocate space for the header structure */
	fitshead = (fitsheader *) calloc(1, sizeof(fitsheader));
	if (!fitshead) {
		error_exit("readfitsheader: failed to allocate space for header\n");
	}

	/* set stream elements */
	fitshead->ipstream = stream;
	fitshead->opstream = stdout;

	/* read the comments into a null terminated, double linked list of fitscomments */
	thecomment = NULL;
	while (!alldone) {
		if (FITS_REC_SIZE != fread(record, sizeof(char), FITS_REC_SIZE, stream)) {
			error_exit("readfitsheader: malformed header\n");
		}
		line = recorddone = 0;
		while (!recorddone) {
			newcomment = (fitscomment *) calloc(1, sizeof(fitscomment));
			if (!newcomment) {
				error_exit("readfitsheader: failed to allocate space for fits header comment\n");
			}
			if (!thecomment) {
				fitshead->basecomment = thecomment = newcomment;
			} else {
				thecomment->next = newcomment;
				newcomment->prev = thecomment;
				thecomment = newcomment;
			}
			strncpy(thecomment->name, record + line * COM_LENGTH, NAME_LENGTH);
			strncpy(thecomment->value, record + line * COM_LENGTH + NAME_LENGTH + 2, VALUE_LENGTH);
			if (!strncmp(thecomment->name, "END     ",  NAME_LENGTH)) {
				alldone = recorddone = 1;
				(newcomment->prev)->next = NULL;
			}
			line++;
			if (line * COM_LENGTH == FITS_REC_SIZE) {
				recorddone = 1;
			}
		}
	}

	/* check/get the required entries */

	/* SIMPLE/XTENSION line */
	bad = 0;
	if (!(thecomment = fitshead->basecomment)) {
		bad = 1;
	}

	if (!bad) {
		if (strncmp(thecomment->name, "SIMPLE  ", NAME_LENGTH)) {
			if (strncmp(thecomment->name, "XTENSION", NAME_LENGTH)) {
				bad = 1;
			} else {
				fitshead->isextension = 1;
			}
		}
	}
	if (!bad) {
		textval = gettextvalue(thecomment);
		if (fitshead->isextension) {
			if (strncmp(textval, "IMAGE", 5)) {
				bad = 1;
			} 
		} else {
			if (strncmp(textval, "T", 1)) {
				bad = 1;
			} 
		}
	}
	if (bad) {
		error_exit("readfitsheader: bad SIMPLE/XTENSION line\n");
	}
	/* now check to see if it is an old imcat file (i.e. no BYTEORDR line) */
	if (strncmp(thecomment->value+34, "imcat", 5)) {
		oldimcatformat = 0;
	} else {
		oldimcatformat = 1;
	}

	/* BITPIX line */
	if (!(thecomment = thecomment->next)) {
		bad = 1;
	}
	if (!bad) {
		if (strcmp(thecomment->name, "BITPIX  ")) {
			bad = 1;
		}
		
	}
	if (!bad) {
		if (1 != sscanf(thecomment->value, "%d", &(fitshead->extpixtype))) {
			bad = 1;
		}
		if (fitshead->extpixtype != UCHAR_PIXTYPE &&
			fitshead->extpixtype != SHORT_PIXTYPE && fitshead->extpixtype != FLOAT_PIXTYPE &&
			fitshead->extpixtype != INT_PIXTYPE && fitshead->extpixtype != DBL_PIXTYPE) {
				bad = 1;
		}
	}
	if (bad) {
		error_exit("readfitsheader: bad BITPIX line\n");
	}

	/* NAXIS line */
	if (!(thecomment = thecomment->next)) {
		bad = 1;
	}
	if (!bad) {
		if (strcmp(thecomment->name, "NAXIS   ")) {
			bad = 1;
		}
		
	}
	if (!bad) {
		if (1 != sscanf(thecomment->value, "%d", &(fitshead->ndim))) {
			bad = 1;
		}
		if (fitshead->ndim >= MAX_FITS_DIM) {
				bad = 1;
		}
	}
	if (bad) {
		error_exit("readfitsheader: bad NAXIS line\n");
	}

	/* axis dimensions */
	for (dim = 0; dim < fitshead->ndim; dim++) {
		if (!(thecomment = thecomment->next)) {
			bad = 1;
		}
		sprintf(string, "NAXIS%-3d", dim + 1);
		if (strcmp(thecomment->name, string)) {
			bad = 1;
		}
		if (1 != sscanf(thecomment->value, "%d", fitshead->n + dim)) {
			bad = 1;
		}
		if (bad) {
			error_exit("readfitsheader: bad NAXIS? line\n");
		}
	}

	/* make fitshead->basecomment point at the next comment */
	fitshead->basecomment = thecomment->next;
	if (fitshead->basecomment) {
		(fitshead->basecomment)->prev = NULL;
	}

	/* deal with extensions */
	if (!(fitshead->isextension)) {
		if (thecomment = getcommentbyname("EXTEND", fitshead)) {
			fitshead->hasextensions = 1;
			removecomment(thecomment, fitshead);
			if (thecomment = getcommentbyname("NEXTEND", fitshead)) {
				fitshead->nextensions = (int) getnumericvalue(thecomment);
				removecomment(thecomment, fitshead);
			}
		}
	} else {
		if (thecomment = getcommentbyname("PCOUNT", fitshead)) {
			fitshead->pcount = (int) getnumericvalue(thecomment);
			if (fitshead->pcount) {
				error_exit("readfitsheader: I can only deal with PCOUNT=0 FITS files\n");
			}
			removecomment(thecomment, fitshead);
		}
		if (thecomment = getcommentbyname("GCOUNT", fitshead)) {
			fitshead->gcount = (int) getnumericvalue(thecomment);
			if (fitshead->gcount != 1) {
				error_exit("readfitsheader: I can only deal with GCOUNT=1 FITS files\n");
			}
			removecomment(thecomment, fitshead);
		}
	}

	/* look for BSCALE, BZERO */
	if (thecomment = getcommentbyname("BSCALE", fitshead)) {
		fitshead->bscale = getnumericvalue(thecomment);
		removecomment(thecomment, fitshead);
		fitshead->bscaling = 1;
	} else {
		fitshead->bscale = 1.0;
	}
	if (thecomment = getcommentbyname("BZERO", fitshead)) {
		fitshead->bzero = getnumericvalue(thecomment);
		removecomment(thecomment, fitshead);
		fitshead->bscaling = 1;
	} else {
		fitshead->bzero = 0.0;
	}

	/* default internal pixtype */
	fitshead->intpixtype = FLOAT_PIXTYPE;

	/* now we set the byte order */
	/* for old imcat format both ip/op controlled solely by IMCATSWAPFITSBYTES */
	if (byteswapenvvar = getenv("IMCATSWAPFITSBYTES")) {
		fitshead->ipbyteorder = NON_NATIVE_BYTE_ORDER;
		fitshead->opbyteorder = NON_NATIVE_BYTE_ORDER;
	} else {
		fitshead->ipbyteorder = NATIVE_BYTE_ORDER;
		fitshead->opbyteorder = NATIVE_BYTE_ORDER;
	}	
	if (!oldimcatformat) {
		/* we will assume it is a NOST conforming file */
		fitshead->ipbyteorder = BIG_ENDIAN_BYTE_ORDER;
		/* unless it explicitly contains the BYTEORDR comment */
		if (thecomment = getcommentbyname("BYTEORDR", fitshead)) {
			if (!strncmp(gettextvalue(thecomment), "LITTLE", 6)) {
				fitshead->ipbyteorder = LITTLE_ENDIAN_BYTE_ORDER;
			}
			removecomment(thecomment, fitshead);
		}
	}

	if (getenv("IMCATCONVERTNANS")) {
		fitshead->convertnans = 1;
	}

	return(fitshead);
}


/*
 * copyfitsheader() makes copy of a header complete with copies of the comment list 
 */
fitsheader *copyfitsheader(fitsheader *srcfits)
{
	fitsheader 	*dstfits;
	fitscomment	*srccom, *dstcom, *newcom;

	dstfits = (fitsheader *) calloc(1, sizeof(fitsheader));
	*dstfits = *srcfits;
	srccom = srcfits->basecomment;
	dstcom = NULL;
	while (srccom) {
		newcom = (fitscomment *) calloc(1, sizeof(fitscomment));
		*newcom = *srccom;
		newcom->next = newcom->prev = NULL;
		if (!dstcom) {
			dstfits->basecomment = dstcom = newcom;
		} else {
			dstcom->next = newcom;
			newcom->prev = dstcom;
			dstcom = newcom;
		}
		srccom = srccom->next;
	}
	return(dstfits);
}


/*
 * writefitsheader() writes header to the stream theheader->opstream
 */
int		writefitsheader(fitsheader *theheader)
{
	char		record[FITS_REC_SIZE], string[COM_LENGTH1];
	int		line, dim;
	fitscomment	*thecomment;
	FILE 		*stream;

	stream = theheader->opstream;

	/* write the required header items */
	line = 0;
	if (theheader->isextension) {
		thecomment = newtextcomment("XTENSION", "'IMAGE   '", "written by IMCAT");
	} else {
		thecomment = newtextcomment("SIMPLE", "T", "written by IMCAT");
	}	
	copycommenttostring(thecomment, record + line++ * COM_LENGTH);
	thecomment = newnumericcomment("BITPIX", (double) theheader->extpixtype, "");	
	copycommenttostring(thecomment, record + line++ * COM_LENGTH);
	thecomment = newnumericcomment("NAXIS", (double) theheader->ndim, "");	
	copycommenttostring(thecomment, record + line++ * COM_LENGTH);
	for (dim = 0; dim < theheader->ndim; dim++) {
		sprintf(string, "NAXIS%d", dim + 1);
		thecomment = newnumericcomment(string, (double) theheader->n[dim], "");	
		copycommenttostring(thecomment, record + line++ * COM_LENGTH);		
	}

	/* EXTENSIONS */
	if (theheader->hasextensions) {
		thecomment = newtextcomment("EXTEND", "T", "file has extensions");
		copycommenttostring(thecomment, record + line++ * COM_LENGTH);
		thecomment = newnumericcomment("NEXTEND", (double) theheader->nextensions, "number of extnesions");
		copycommenttostring(thecomment, record + line++ * COM_LENGTH);
	}
	if (theheader->isextension) {
		thecomment = newnumericcomment("PCOUNT", 0.0, "no random parameters");
		copycommenttostring(thecomment, record + line++ * COM_LENGTH);
		thecomment = newnumericcomment("GCOUNT", 1.0, "single group");
		copycommenttostring(thecomment, record + line++ * COM_LENGTH);
	}

	/* BSCALE and BZERO */
	if (theheader->extpixtype == FLOAT_PIXTYPE || theheader->extpixtype == DBL_PIXTYPE) {
		theheader->bscaling = 0;
	}
	if (theheader->bscaling) {
		thecomment = newnumericcomment("BSCALE", theheader->bscale, "");	
		copycommenttostring(thecomment, record + line++ * COM_LENGTH);
		thecomment = newnumericcomment("BZERO", theheader->bzero, "");	
		copycommenttostring(thecomment, record + line++ * COM_LENGTH);
	}

	/* BYTE ORDER */
	if (theheader->opbyteorder == BIG_ENDIAN_BYTE_ORDER) {
		thecomment = newtextcomment("BYTEORDR", "BIG_ENDIAN", "SunOS, solaris etc. byte order");
	} else {
		thecomment = newtextcomment("BYTEORDR", "LITTLE_ENDIAN", "OSF1, Linux etc byte order");
	}
	copycommenttostring(thecomment, record + line++ * COM_LENGTH);
	
	/* output the linked list */
	thecomment = theheader->basecomment;
	while (thecomment) {
		copycommenttostring(thecomment, record + line++ * COM_LENGTH);		
		thecomment = thecomment->next;
		if (line * COM_LENGTH == FITS_REC_SIZE) {
			if (FITS_REC_SIZE != fwrite(record, sizeof(char), FITS_REC_SIZE, stream)) {
				error_exit("writefitsheader: fwrite failed\n");
			}
			line = 0;
		}
	}

	/* output the END comment */
	sprintf(string, "%-80.80s", "END");
	strncpy(record + line++ * COM_LENGTH, string, COM_LENGTH);
	while (line * COM_LENGTH < FITS_REC_SIZE) {
		sprintf(string, "%-80.80s", "");
		strncpy(record + line++ * COM_LENGTH, string, COM_LENGTH);
	}
	if (FITS_REC_SIZE != fwrite(record, sizeof(char), FITS_REC_SIZE, stream)) {
		error_exit("writefitsheader: fwrite failed\n");
	}
	
	return (1);
}


/*
 * writefitstail() writes tail of last 2880 byte disk record
 */
int		writefitstail(fitsheader *theheader)
{
	int 	i, ic, pad;
	char	record[2880];
	long	length;
	FILE 	*stream;

	stream = theheader->opstream;

	switch (theheader->extpixtype) {
		case UCHAR_PIXTYPE:
			length = (long) sizeof(unsigned char);
			break;
		case SHORT_PIXTYPE:
			length = (long) sizeof(short);
			break;
		case FLOAT_PIXTYPE:
			length = (long) sizeof(float);
			break;
		case INT_PIXTYPE:
			length = (long) sizeof(int);
			break;
		case DBL_PIXTYPE:
			length = (long) sizeof(double);
			break;
		default:
			error_exit("write_fits_tail: illegal pixtype\n");
			break;
	}
	for (i = 0; i < theheader->ndim; i++) {
		length *= (long) (theheader->n)[i];
	}
	pad = 2880 - length + 2880 * (length / 2880);
	for (ic = 0; ic < pad; ic++) {
		record[ic] = ' ';
	}
	return(fwrite(record, sizeof(char), pad, stream));	
}


int	 	copycommenttostring(fitscomment *thecomment, char *string)
{
	char	tempstring[COM_LENGTH1];

	sprintf(tempstring, "%-8.8s= %-70.70s", thecomment->name, thecomment->value);
	strncpy(string, tempstring, COM_LENGTH);
}



fitscomment	*newtextcomment(char *name, char *value, char *comment)
{
	fitscomment	*thecomment;


	if (strlen(name) > NAME_LENGTH) {
		error_exit("newtextcomment: name too long\n");
	}
	if (!(thecomment = (fitscomment *) calloc(1, sizeof(fitscomment)))) {
		error_exit("newtextcomment: memory allocation failure\n");
	}
	sprintf(thecomment->name, "%-8.8s", name);
	if (comment) {
		sprintf(thecomment->value, "%20s / %-47.47s", value, comment);
	} else {
		sprintf(thecomment->value, "%-70.70s", value);
	}
	return(thecomment);
}


fitscomment	*newnumericcomment(char *name, double value, char *comment)
{
	fitscomment	*thecomment;


	if (strlen(name) > NAME_LENGTH) {
		error_exit("newnumericcomment: name too long\n");
	}
	if (!(thecomment = (fitscomment *) calloc(1, sizeof(fitscomment)))) {
		error_exit("newnumericcomment: memory allocation failure\n");
	}
	sprintf(thecomment->name, "%-8.8s", name);
	if (comment) {
		sprintf(thecomment->value, "%20.8lg / %-47.47s", value, comment);
	} else {
		sprintf(thecomment->value, "%20.8lg", value);
	}
	return(thecomment);
}



fitscomment	*getcommentbyname(char *name, fitsheader *theheader)
{
	fitscomment	*thecomment;
	char		string[NAME_LENGTH1];

	if (strlen(name) > NAME_LENGTH) {
		error_exit("getcommentbyname: name too long\n");
	}
	sprintf(string, "%-8.8s", name);
	thecomment = theheader->basecomment;
	while (thecomment) {
		if (!strncmp(thecomment->name, string, NAME_LENGTH)) {
			return(thecomment);
		}
		thecomment = thecomment->next;
	}
	return(NULL);
}


double		getnumericvalue(fitscomment *thecomment)
{
	double	val;

	if (1 != sscanf(thecomment->value, "%lg", &val)) {
		error_exit("getnumericvalue: can't decipher value\n");
	}
	return(val);
}


char		*gettextvalue(fitscomment *thecomment)
{
	char	*val, *tmpstr;
	int	pos = 0;
	
	tmpstr = (char *) calloc(VALUE_LENGTH1, sizeof(char));
	strncpy(tmpstr, thecomment->value, VALUE_LENGTH);
	while (!strncmp(tmpstr, " ", 1)) {
		tmpstr++;
	}
	if (strncmp(tmpstr, "'", 1)) {
		val = strtok(tmpstr, " ");
	} else {
		val = strtok(tmpstr, "'");
	}
	return(val);
}


void		readfitsline(void *f, fitsheader *theheader)
{	
	int	N1;
	FILE	*theipf;
	
	N1 = theheader->n[0];
	theipf = theheader->ipstream;

	if (theheader->extpixtype == theheader->intpixtype && theheader->ipbyteorder == NATIVE_BYTE_ORDER &&
			!(theheader->convertnans && (theheader->extpixtype == FLOAT_PIXTYPE || theheader->extpixtype == DBL_PIXTYPE))) {
		fread(f, pixsize(theheader->extpixtype), N1, theipf);
	} else {
		updatelinebuffer(theheader);
		fread(theheader->linebuffer, pixsize(theheader->extpixtype), N1, theipf);
		if (theheader->ipbyteorder == NON_NATIVE_BYTE_ORDER) {
			byteswapline((void *) (theheader->linebuffer), N1, pixsize(theheader->extpixtype));
		}
		if (theheader->convertnans && (theheader->extpixtype == FLOAT_PIXTYPE || theheader->extpixtype == DBL_PIXTYPE)) {
			convertnanstomagic((void *) (theheader->linebuffer), N1, pixsize(theheader->extpixtype));
		}
		convertarray((char *) theheader->linebuffer, (char *) f, theheader->extpixtype, theheader->intpixtype, N1, 
				theheader->bscaling, theheader->bscale, theheader->bzero);
	}
}


void		writefitsline(void *f, fitsheader *theheader)
{	
	int	N1;
	FILE	*theopf;
	double	bscale, bzero;
	
	N1 = theheader->n[0];
	theopf = theheader->opstream;
	if (theheader->bscaling) {
		bscale = 1.0 / theheader->bscale;
		bzero = - theheader->bzero / theheader->bscale;
	}

	if (theheader->extpixtype == theheader->intpixtype && theheader->opbyteorder == NATIVE_BYTE_ORDER && 
			!(theheader->convertnans && (theheader->extpixtype == FLOAT_PIXTYPE || theheader->extpixtype == DBL_PIXTYPE))) {
		fwrite(f, pixsize(theheader->extpixtype), N1, theopf);
	} else {
		updatelinebuffer(theheader);
		convertarray((char *) f, (char *) theheader->linebuffer, theheader->intpixtype, theheader->extpixtype, N1, 
				theheader->bscaling, bscale, bzero);
		if (theheader->convertnans && (theheader->extpixtype == FLOAT_PIXTYPE || theheader->extpixtype == DBL_PIXTYPE)) {
			convertmagictonans((void *) (theheader->linebuffer), N1, pixsize(theheader->extpixtype));
		}
		if (theheader->opbyteorder == NON_NATIVE_BYTE_ORDER) {
			byteswapline((void *) (theheader->linebuffer), N1, pixsize(theheader->extpixtype));
		}
		fwrite(theheader->linebuffer, pixsize(theheader->extpixtype), N1, theopf);
	}
}

void		readfitsplane(void **f, fitsheader *theheader)
{
	int	i;

	for (i = 0; i < theheader->n[1]; i++) {
		readfitsline(f[i], theheader);
	}
}

void		writefitsplane(void **f, fitsheader *theheader)
{
	int	i;

	for (i = 0; i < theheader->n[1]; i++) {
		writefitsline(f[i], theheader);
	}
}

void		readfitscube(void ***f, fitsheader *theheader)
{
	int	i;

	for (i = 0; i < theheader->n[2]; i++) {
		readfitsplane(f[i], theheader);
	}
}

void		writefitscube(void ***f, fitsheader *theheader)
{
	int	i;

	for (i = 0; i < theheader->n[2]; i++) {
		writefitsplane(f[i], theheader);
	}
}

void		removecomment(fitscomment *thecomment, fitsheader *fitshead)
{
	fitscomment *com;

	if (fitshead->basecomment == thecomment) {
		fitshead->basecomment = thecomment->next;
		return;
	}
	com = fitshead->basecomment;
	while (com) {
		if (com == thecomment) {
			if (thecomment->prev) {
				(thecomment->prev)->next = thecomment->next;
			}
			if (thecomment->next) {
				(thecomment->next)->prev = thecomment->prev;
			}
			com = NULL;
			return;
		}
		com = com->next;
	}
}


void		removenamedcomments(char *name, fitsheader *fitshead)
{
	fitscomment *com;	

	while (com = getcommentbyname(name, fitshead)) {
		removecomment(com, fitshead);
	}
}


int		pixsize(int pixtype)
{
	switch (pixtype) {
		case UCHAR_PIXTYPE:
			return(sizeof(unsigned char));
			break; 
		case SHORT_PIXTYPE:
			return(sizeof(short));
			break; 
		case INT_PIXTYPE:
			return(sizeof(int));
			break; 
		case FLOAT_PIXTYPE:
			return(sizeof(float));
			break; 
		case DBL_PIXTYPE:
			return(sizeof(double));
			break;
		default:
			error_exit("pixsize: bad pixtype\n");
	}
}


/* higher level routines */

int		read2Dfloatimage(float ***f, int *N1, int *N2, fitsheader **fits, FILE *stream)
{
	int	y;

	*fits = readfitsheader(stream);
	if ((*fits)->ndim != 2) {
		error_exit("read2Dfloatimage: image not 2 dimensional\n");
	}
	*N1 = (*fits)->n[0];
	*N2 = (*fits)->n[1];
	allocFloatArray(f, *N1, *N2);
	for (y = 0; y < *N2; y++) {
		readfitsline((void *) (*f)[y], *fits);
	}
}

int		read2Dfloatimage_shm(float ***f, int *N1, int *N2, fitsheader **fits, FILE *stream)
{
	error_exit("read2Dfloatimage_shm: not implemented right now\n");
}

int		write2Dfloatimage(float **f, fitsheader *fits)
{
	int	y;

	if (fits->ndim != 2) {
		error_exit("write2Dfloatimage: image not 2 dimensional\n");
	}
	writefitsheader(fits);
	for (y = 0; y < fits->n[1]; y++) {
		writefitsline((void *) f[y], fits);
	}
	writefitstail(fits);		
}


fitsheader 	*newfitsheader(int ndim, int *dim, int extpixtype)
{
	fitsheader 	*fits;
	char		*byteswapenvvar;
	int		idim;

	fits = (fitsheader *) calloc(1, sizeof(fitsheader));
	fits->extpixtype = extpixtype;
	fits->intpixtype = FLOAT_PIXTYPE;
	fits->ndim = ndim;
	for (idim = 0; idim < ndim; idim++) {
		fits->n[idim] = dim[idim];
	}
	fits->bscaling = 0;
	fits->ipstream = stdin;
	fits->opstream = stdout;
	if (byteswapenvvar = getenv("IMCATSWAPFITSBYTES")) {
		fits->opbyteorder = NON_NATIVE_BYTE_ORDER;
	} else {
		fits->opbyteorder = NATIVE_BYTE_ORDER;
	}
	if (getenv("IMCATCONVERTNANS")) {
		fits->convertnans = 1;
	}
	return(fits);
}

fitsheader 	*new2Dfitsheader(int N1, int N2, int extpixtype)
{
	int	dim[2];

	dim[0] = N1;
	dim[1] = N2;
	return(newfitsheader(2, dim, extpixtype));
}


void	argsToString(int argc, char **argv, char *string)
{
	int	i, stringlen, arglen;

	if (strlen(argv[0]) > COM_LENGTH) {
		strncpy(string, argv[0], COM_LENGTH);
		string[COM_LENGTH] = '\0';
	} else {
		strcpy(string, argv[0]);
	}
	for (i = 1; i < argc; i++) {
		stringlen = strlen(string);
		arglen = strlen(argv[i]);
		if (stringlen + 1 + arglen <= COM_LENGTH) {
			strcat(string, " ");
			strcat(string, argv[i]);
		} else {
			break;
		}
	}
}

void	add_comment(int argc, char **argv, fitsheader *fitshead)
{
	char		argstring[COM_LENGTH1];
	fitscomment	*thecomment;

	argsToString(argc, argv, argstring);
	thecomment = newtextcomment("HISTORY", argstring, NULL);
	appendcomment(thecomment, fitshead);
}


void	appendcomment(fitscomment *newcomment, fitsheader *fitshead)
{
	fitscomment	*thecomment;

	if (!fitshead->basecomment) {
		fitshead->basecomment = newcomment;	
		return;
	}
	thecomment = fitshead->basecomment;
	while (thecomment->next) {
		thecomment = thecomment->next;
	}
	thecomment->next = newcomment;
}


void	prependcomment(fitscomment *newcomment, fitsheader *fitshead)
{
	if (fitshead->basecomment) {
		newcomment->next = fitshead->basecomment;
		(fitshead->basecomment)->prev = newcomment;
	}
	fitshead->basecomment = newcomment;
}


void		setextpixtype(fitsheader *fits, int pixtype)
{
	fits->extpixtype = pixtype;
	fits->linebuffer = NULL;
}


void		set2Dimagesize(fitsheader *fits, int N1, int N2)
{
	fits->ndim = 2;
	fits->n[0] = N1;
	fits->n[1] = N2;
	fits->linebuffer = NULL;
}


int		skiplines(fitsheader *fits, int nlines)
{
	long	offset;
	struct	stat	st;
	
	/* first we see if the file is regular */
	fstat(fileno(fits->ipstream), &st);
	if ((st.st_mode & S_IFMT) == S_IFREG) {
		offset = (long) (fits->n[0] * pixsize(fits->extpixtype) * nlines);
		return(fseek(fits->ipstream, offset, SEEK_CUR));
	} else {
		updatelinebuffer(fits);
		while (nlines--) {
			fread((void *) fits->linebuffer, pixsize(fits->extpixtype), fits->n[0], fits->ipstream); 
		}
	}
}

int		updatelinebuffer(fitsheader *theheader) 
{
	int	N1;

	N1 = theheader->n[0];
	if (!theheader->linebuffer || theheader->linebuffersize != N1 * pixsize(theheader->extpixtype)) {
		if (theheader->linebuffer) {
			free(theheader->linebuffer);
		}
		theheader->linebuffer = (char *) calloc(N1 * pixsize(theheader->extpixtype), sizeof(char));
		theheader->linebuffersize = N1 * pixsize(theheader->extpixtype);
	}
}


int	byteswapline(void *data, int nel, int pixsize)
{
	static char	tmpchar[8];
	int		i, b;
	char		*caddr;

        caddr = (char *) data;
        for (i = 0; i < nel; i++) {
                for (b = 0 ; b < pixsize; b++) {
                        tmpchar[b] = caddr[b];
                }
                for (b = 0 ; b < pixsize; b++) {
                        caddr[b] = tmpchar[pixsize - b - 1];
                }
                caddr += pixsize;
        }

}


int	convertmagictonans(void *data, int nel, int pixsize)
{
	_Dconst _Nan = {{INIT(_DNAN)}};
	double	dnan = _Nan._D;
	float	fnan = (float) dnan;
	int	i;

	switch (pixsize) {
		case 4:
			for (i = 0; i < nel; i++) {
				if (((float *) data)[i] == FLOAT_MAGIC) {
					((float *) data)[i] = fnan;
				}
			}
			break;
		case 8:
			for (i = 0; i < nel; i++) {
				if (((double *) data)[i] == DBL_MAGIC) {
					((double *) data)[i] = dnan;
				}
			}
			break;
		default:
			error_exit("convertmagictonans: illegal pixsize\n");
	}		
}

int	convertnanstomagic(void *data, int nel, int pixsize)
{
	_Dconst _Nan = {{INIT(_DNAN)}};
	double	dnan = _Nan._D;
	float	fnan = (float) dnan;
	int	i;

	switch (pixsize) {
		case 4:
			for (i = 0; i < nel; i++) {
				if (isnan((double) (((float *) data)[i]))) {
					((float *) data)[i] = FLOAT_MAGIC;
				}
			}
			break;
		case 8:
			for (i = 0; i < nel; i++) {
				if (isnan(((double *) data)[i])) {
					((double *) data)[i] = DBL_MAGIC;
				}
			}
			break;
		default:
			error_exit("convertnanstomagic: illegal pixsize\n");
	}		
}
