/*** File libwcs/hput.c
 *** July 17, 1997
 *** By Doug Mink

 * Module:	hput.c (Put FITS Header parameter values)
 * Purpose:	Implant values for parameters into FITS header string
 * Subroutine:	hputi2 (hstring,keyword,ival) sets integer*2 ival
 * Subroutine:	hputi4 (hstring,keyword,ival) sets int ival
 * Subroutine:	hputr4 (hstring,keyword,rval) sets real*4 rval
 * Subroutine:	hputr8 (hstring,keyword,dval) sets real*8 dval
 * Subroutine:	hputd8 (hstring,keyword,ndec,dval) sets real*8 dval
 * Subroutine:	hputra (hstring,keyword,lval) sets right ascension as string
 * Subroutine:	hputdec (hstring,keyword,lval) sets declination as string
 * Subroutine:	hputl  (hstring,keyword,lval) sets logical lval
 * Subroutine:	hputs  (hstring,keyword,cval) sets character string adding ''
 * Subroutine:	hputc  (hstring,keyword,cval) sets character string cval
 * Subroutine:	hdel   (hstring,keyword) deletes entry for keyword keyword
 * Subroutine:	hchange (hstring,keyword1,keyword2) changes keyword for entry
 * Subroutine:	hputcom (hstring,keyword,comment) sets comment for parameter keyword
 * Subroutine:	ra2str (out, ra, ndec) converts RA from degrees to string
 * Subroutine:	dec2str (out, dec, ndec) converts Dec from degrees to string
 * Subroutine:	deg2str (out, deg, ndec) converts degrees to string

 * Copyright:   1995, 1996 Smithsonian Astrophysical Observatory
 *              You may do anything you like with this file except remove
 *              this copyright.  The Smithsonian Astrophysical Observatory
 *              makes no representations about the suitability of this
 *              software for any purpose.  It is provided "as is" without
 *              express or implied warranty.
 */
#include <string.h>             /* NULL, strlen, strstr, strcpy */
#include <stdio.h>
#include <stdlib.h>
#include "fitshead.h"

static int verbose=0;	/* Set to 1 to print error messages and other info */
void hputc();


/*  HPUTI4 - Set int keyword = ival in FITS header string */

void
hputi4 (hstring,keyword,ival)

  char *hstring;	/* character string containing FITS-style header
			   information in the format
			   <keyword>= <value> {/ <comment>}
			   each entry is padded with spaces to 80 characters */

  char *keyword;		/* character string containing the name of the variable
			   to be returned.  hput searches for a line beginning
			   with this string, and if there isn't one, creates one.
		   	   The first 8 characters of keyword must be unique. */
  int ival;		/* int number */
{
    char value[30];

    /* Translate value from binary to ASCII */
    sprintf (value,"%d",ival);

    /* Put value into header string */
    hputc (hstring,keyword,value);

    /* Return to calling program */
    return;
}


/*  HPUTI2 - Set short keyword = ival in FITS header string */

void
hputi2 (hstring,keyword,ival)

  char *hstring;	/* FITS header string */
  char *keyword;		/* Keyword name */
  short ival;		/* short number */

{
    char value[30];

    /* Translate value from binary to ASCII */
    sprintf (value,"%d",ival);

    /* Put value into header string */
    hputc (hstring,keyword,value);

    /* Return to calling program */
    return;
}


/*  HPUTR4 - Set float keyword = rval in FITS header string */

void
hputr4 (hstring,keyword,rval)

char *hstring;		/* FITS header string */
char *keyword;		/* Keyword name */
float rval;		/* float number */
{
    char value[30];

    /* Translate value from binary to ASCII */
    sprintf (value,"%f",rval);

    /* Put value into header string */
    hputc (hstring,keyword,value);

    /* Return to calling program */
    return;
}


/*  HPUTR8 - Set double keyword = dval in FITS header string */

void
hputr8 (hstring,keyword,dval)

char	*hstring;	/* FITS header string */
char	*keyword;	/* Keyword name */
double	dval;		/* double number */
{
    char value[30];

    /* Translate value from binary to ASCII */
    sprintf (value,"%g",dval);

    /* Put value into header string */
    hputc (hstring,keyword,value);

    /* Return to calling program */
    return;
}


/*  HPUTNR8 - Set double keyword = dval in FITS header string */

void
hputnr8 (hstring,keyword,ndec,dval)

char	*hstring;	/* FITS header string */
char	*keyword;	/* Keyword name */
int	ndec;		/* Number of decimal places to print */
double	dval;		/* double number */
{
    char value[30];
    char format[8];

    /* Translate value from binary to ASCII */
    if (ndec < 0)
	sprintf (value,"%g",dval);
    else {
	sprintf (format,"%%.%df",ndec);
	sprintf (value,format,dval);
	}

    /* Put value into header string */
    hputc (hstring,keyword,value);

    /* Return to calling program */
    return;
}


/*  HPUTRA - Set double keyword = hh:mm:ss.sss in FITS header string */

void
hputra (hstring,keyword, ra)

char *hstring;		/* FITS header string */
char *keyword;		/* Keyword name */
double ra;		/* Right ascension in degrees */
{
    char value[30];

    /* Translate value from binary to ASCII */
    ra2str (value, ra, 3);

    /* Put value into header string */
    hputs (hstring,keyword,value);

    /* Return to calling program */
    return;
}


/*  HPUTDEC - Set double keyword = dd:mm:ss.sss in FITS header string */

void
hputdec (hstring, keyword, dec)

char *hstring;		/* FITS header string */
char *keyword;		/* Keyword name */
double dec;		/* Declination in degrees */
{
    char value[30];

    /* Translate value from binary to ASCII */
    dec2str (value, dec, 2);

    /* Put value into header string */
    hputs (hstring,keyword,value);

    /* Return to calling program */
   return;
}



/*  HPUTL - Set keyword = F if lval=0, else T, in FITS header string */

void
hputl (hstring, keyword,lval)

char *hstring;		/* FITS header */
char *keyword;		/* Keyword name */
int lval;		/* logical variable (0=false, else true) */
{
    char value[8];

    /* Translate value from binary to ASCII */
    if (lval)
	strcpy (value, "T");
    else
	strcpy (value, "F");

    /* Put value into header string */
    hputc (hstring,keyword,value);

    /* Return to calling program */
    return;
}


/*  HPUTS - Set character string keyword = 'cval' in FITS header string */

void
hputs (hstring,keyword,cval)

char *hstring;	/* FITS header */
char *keyword;	/* Keyword name */
char *cval;	/* character string containing the value for variable
		   keyword.  trailing and leading blanks are removed.  */
{
    char squot = 39;
    char value[70];
    int lcval;

    /*  find length of variable string */

    lcval = strlen (cval);
    if (lcval > 67)
	lcval = 67;

    /* Put quotes around string */
    value[0] = squot;
    strncpy (&value[1],cval,lcval);
    value[lcval+1] = squot;
    value[lcval+2] = 0;

    /* Put value into header string */
    hputc (hstring,keyword,value);

    /* Return to calling program */
    return;
}


/*  HPUTC - Set character string keyword = value in FITS header string */

void
hputc (hstring,keyword,value)

char *hstring;
char *keyword;
char *value;	/* character string containing the value for variable
		   keyword.  trailing and leading blanks are removed.  */
{
	char squot = 39;
	char line[100];
	char newcom[50];
	char *v, *vp, *v1, *v2, *q1, *q2, *c1, *ve;
	char *ksearch();
	int lkeyword, lcom, lval, lc;

/*  find length of keyword and value */
	lkeyword = strlen (keyword);
	lval = strlen (value);

/*  If COMMENT or HISTORY, always add it just before the END */
	if (lkeyword == 7 && (strncmp (keyword,"COMMENT",7) == 0 ||
	    strncmp (keyword,"HISTORY",7) == 0)) {

	/* Find end of header */
	    v1 = ksearch (hstring,"END");
	    v2 = v1 + 80;

	/* Move END down one line */
	    strncpy (v2, v1, 80);

	/* Insert keyword */
	    strncpy (v1,keyword,7);

	/* Pad with spaces */
	    for (vp = v1+lkeyword; vp < v2; vp++)
		*vp = ' ';

	/* Insert comment */
	    strncpy (v1+9,value,lval);
	    return;
	    }

/* Otherwise search for keyword */
	else
	    v1 = ksearch (hstring,keyword);

/*  Find end of header */
	ve = ksearch (hstring,"END");

/*  If parameter is not found, create a space for it */
	if (v1 == NULL) {
	    v1 = ve;
	    v2 = v1 + 80;
	    strncpy (v2, ve, 80);
	    lcom = 0;
	    newcom[0] = 0;
	    }

/*  Otherwise, extract the entry for this keyword from the header */
	else {
	    strncpy (line, v1, 80);
	    line[80] = 0;
	    v2 = v1 + 80;

/*  check for quoted value */
	    q1 = strchr (line, squot);
	    if (q1 != NULL)
		q2 = strchr (q1+1,squot);
	    else
		q2 = line;

/*  extract comment and remove trailing spaces */

	    c1 = strchr (q2,'/');
	    if (c1 != NULL) {
		lcom = 80 - (c1 -line);
		strncpy (newcom, c1+1, lcom);
		vp = newcom + lcom - 1;
		while (vp-- > newcom && *vp == ' ')
		    *vp = 0;
		lcom = strlen (newcom);
		}
	    else {
		newcom[0] = 0;
		lcom = 0;
		}
	    }

/* Fill new entry with spaces */
	for (vp = v1; vp < v2; vp++)
	    *vp = ' ';

/*  Copy keyword to new entry */
	strncpy (v1, keyword, lkeyword);

/*  Add parameter value in the appropriate place */
	vp = v1 + 8;
	*vp = '=';
	vp = v1 + 9;
	*vp = ' ';
	vp = vp + 1;
	if (*value == squot) {
	    strncpy (vp, value, lval);
	    if (lval+12 > 31)
		lc = lval + 12;
	    else
		lc = 30;
	    }
	else {
	    vp = v1 + 30 - lval;
	    strncpy (vp, value, lval);
	    lc = 30;
	    }

/* Add comment in the appropriate place */
	if (lcom > 0) {
	    if (lc+2+lcom > 80)
		lcom = 78 - lc;
	    vp = v1 + lc + 1;     /* Jul 16 1997: was vp = v1 + lc * 2 */
	    strncpy (vp, "/", 1);
	    vp = vp + 1;
	    strncpy (vp, newcom, lcom);
	    for (v = vp + lcom; v < v2; v++)
		*v = ' ';
	    }

	if (verbose) {
	    if (lcom > 0)
		printf ("HPUT: %s  = %s  / %s\n",keyword, value, newcom);
	    else
		printf ("HPUT: %s  = %s\n",keyword, value);
	    }

	return;
}


/*  HPUTCOM - Set comment for keyword or on line in FITS header string */

void
hputcom (hstring,keyword,comment)

  char *hstring;
  char *keyword;
  char *comment;
{
	char squot;
	char line[100];
	int lkeyword, lcom;
	char *vp, *v1, *v2, *c0, *c1, *q1, *q2;
	char *ksearch();

	squot = 39;

/*  Find length of variable name */
	lkeyword = strlen (keyword);

/*  If COMMENT or HISTORY, always add it just before the END */
	if (lkeyword == 7 && (strncmp (keyword,"COMMENT",7) == 0 ||
	    strncmp (keyword,"HISTORY",7) == 0)) {

	/* Find end of header */
	    v1 = ksearch (hstring,"END");
	    v2 = v1 + 80;
	    strncpy (v2, v1, 80);

	/*  blank out new line and insert keyword */
	    for (vp = v1; vp < v2; vp++)
		*vp = ' ';
	    strncpy (v1, keyword, lkeyword);
	    }

/* search header string for variable name */
	else {
	    v1 = ksearch (hstring,keyword);
	    v2 = v1 + 80;

	/* if parameter is not found, return without doing anything */
	    if (v1 == NULL) {
		if (verbose)
		    printf ("HPUTCOM: %s not found\n",keyword);
		return;
		}

	/* otherwise, extract entry for this variable from the header */
	    strncpy (line, v1, 80);

	/* check for quoted value */
	    q1 = strchr (line,squot);
	    if (q1 != NULL)
		q2 = strchr (q1+1,squot);
	    else
		q2 = NULL;

	    if (q2-line < 31)
		c0 = v1 + 31;
	    else
		c0 = q2 + 2;

	    strncpy (c0, "/ ",2);
	    }

/* create new entry */
	lcom = strlen (comment);

	if (lcom > 0) {
	    c1 = c0 + 2;
	    if (c1+lcom > v2)
		lcom = v2 - c1;
	    strncpy (c1, comment, lcom);
	    }

	if (verbose) {
	    printf ("HPUTCOM: %s / %s\n",keyword,comment);
	    }
}


/*  HDEL - Set character string keyword = value in FITS header string
 *	    returns 1 if entry deleted, else 0
 */

int
hdel (hstring,keyword)

char *hstring;		/* FITS header */
char *keyword;		/* Keyword of entry to be deleted */
{
    char *v, *v1, *v2, *ve;
    char *ksearch();

    /* Search for keyword */
    v1 = ksearch (hstring,keyword);

    /*  If keyword is not found, return header unchanged */
    if (v1 == NULL) {
	return (0);
	}

    /*  Find end of header */
    ve = ksearch (hstring,"END");

    /* Shift rest of header up one line */
    for (v = v1; v < ve; v = v + 80) {
	v2 = v + 80;
	strncpy (v, v2, 80);
	}

    /* Cover former last line with spaces */
    v2 = ve + 80;
    for (v = ve; v < v2; v++)
	*v = ' ';

    return (1);
}


/*  HCHANGE - Changes keyword for entry from keyword1 to keyword2 in FITS
              header string
 *	      returns 1 if entry changed, else 0
 */

int
hchange (hstring, keyword1, keyword2)

char *hstring;		/* FITS header */
char *keyword1;		/* Keyword to be changed */
char *keyword2;		/* New keyword name */
{
    char *v, *v1, *v2;
    int lv2, i;
    char *ksearch();

    /* Search for keyword */
    v1 = ksearch (hstring,keyword1);

    /*  If keyword is not found, return header unchanged */
    if (!v1)
	return (0);

    else {
	lv2 = strlen (keyword2);
	v = v1;
	v2 = keyword2;
	for (i = 0; i < 8; i++) {
	    if (i < lv2)
		v[i] = v2[i];
	    else
		v[i] = ' ';
	    }
	}

    return (1);
}


/* sprint the right ascension ra in sexagesimal format into out[] */

void
ra2str (string, ra, ndec)

char	*string;	/* Character string (returned) */
double	ra;		/* Right ascension in degrees */
int	ndec;		/* Number of decimal places in seconds */

{
    double a,b;
    double seconds;
    int hours;
    int minutes;
    int isec;

    a = ra / 15.0;

    /* Convert to hours */
    hours = (int) a;

    /* Compute minutes */
    b =  (a - (double)hours) * 60.0;
    minutes = (int) b;

    /* Compute seconds */
    seconds = (b - (double)minutes) * 60.0;
    isec = (int)(seconds + 0.5);

    if (ndec > 5) {
	if (seconds > 59.999999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	if (hours > 23)
	    hours = hours - 24;
	(void) sprintf (string,"%2d:%02d:%09.6f",hours,minutes,seconds);
	}
    else if (ndec > 4) {
	if (seconds > 59.99999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	if (hours > 23)
	    hours = hours - 24;
	(void) sprintf (string,"%2d:%02d:%08.5f",hours,minutes,seconds);
	}
    else if (ndec > 3) {
	if (seconds > 59.9999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	if (hours > 23)
	    hours = hours - 24;
	(void) sprintf (string,"%2d:%02d:%07.4f",hours,minutes,seconds);
	}
    else if (ndec > 2) {
	if (seconds > 59.999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	if (hours > 23)
	    hours = hours - 24;
	(void) sprintf (string,"%2d:%02d:%06.3f",hours,minutes,seconds);
	}
    else if (ndec > 1) {
	if (seconds > 59.99) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	if (hours > 23)
	    hours = hours - 24;
	(void) sprintf (string,"%2d:%02d:%05.2f",hours,minutes,seconds);
	}
    else if (ndec > 0) {
	if (seconds > 59.9) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	if (hours > 23)
	    hours = hours - 24;
	(void) sprintf (string,"%2d:%02d:%04.1f",hours,minutes,seconds);
	}
    else if (ndec > -1) {
	if (isec > 59) {
	    isec = 0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	if (hours > 23)
	    hours = hours - 24;
	(void) sprintf (string,"%2d:%02d:%04.1f",hours,minutes,seconds);
	}
    return;
}


/* sprint the variable a in sexagesimal format into string[]. */

void
dec2str (string, dec, ndec)

char	*string;	/* Character string (returned) */
double	dec;		/* Declination in degrees */
int	ndec;		/* Number of decimal places in arcseconds */

{
    double a,b;
    double seconds;
    char sign;
    int degrees;
    int minutes;
    int isec;

    a = dec;

    /* Set sign and do all the rest with a positive */
    if (a < 0) {
	sign = '-';
	a = -a;
	}
    else
	sign = '+';

    /* Convert to degrees */
    degrees = (int) a;

    /* Compute minutes */
    b =  (a - (double)degrees) * 60.0;
    minutes = (int) b;

    /* Compute seconds */
    seconds = (b - (double)minutes) * 60.0;
    isec = (int)(seconds + 0.5);

    if (ndec > 5) {
	if (seconds > 59.999999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (string,"%c%02d:%02d:%09.6f",sign,degrees,minutes,seconds);
	}
    else if (ndec > 4) {
	if (seconds > 59.99999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (string,"%c%02d:%02d:%08.5f",sign,degrees,minutes,seconds);
	}
    else if (ndec > 3) {
	if (seconds > 59.9999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (string,"%c%02d:%02d:%07.4f",sign,degrees,minutes,seconds);
	}
    else if (ndec > 2) {
	if (seconds > 59.999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (string,"%c%02d:%02d:%06.3f",sign,degrees,minutes,seconds);
	}
    else if (ndec > 1) {
	if (seconds > 59.99) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (string,"%c%02d:%02d:%05.2f",sign,degrees,minutes,seconds);
	}
    else if (ndec > 0) {
	if (seconds > 59.9) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (string,"%c%02d:%02d:%04.1f",sign,degrees,minutes,seconds);
	}
    else if (ndec > -1) {
	if (isec > 59) {
	    isec = 0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (string,"%c%02d:%02d:%04.1f",sign,degrees,minutes,seconds);
	}
   return;
}


/* sprint the variable a in sexagesimal format into string[]. */

void
deg2str (string, deg, ndec)

char	*string;	/* Character string (returned) */
double	deg;		/* Angle in degrees */
int	ndec;		/* Number of decimal places in degree string */

{
    char degform[8];
    int field;
    double deg1;

    field = ndec + 4;
    if (ndec > 0) {
	sprintf (degform, "%%%d.%df", field, ndec);
	sprintf (string, degform, deg);
	}
    else {
	sprintf (degform, "%%%4d", field);
	sprintf (string, degform, (int)deg);
	}
    deg1 = atof (string);
    if (deg1 >= 360.0)
	deg1 = deg1 - 360.0;
    else if (deg1 <= -180.0)
	deg1 = deg1 + 360.0;
    if (ndec > 0)
	sprintf (string, degform, deg);
    else
	sprintf (string, degform, (int)deg);
    return;
}

/* Dec 14 1995	Original subroutines

 * Feb  5 1996	Added HDEL to delete keyword entry from FITS header
 * Feb  7 1996	Add EOS to LINE in HPUTC
 * Feb 21 1996	Add RA2STR and DEC2STR string routines
 * Jul 19 1996	Add HPUTRA and HPUTDEC
 * Jul 22 1996	Add HCHANGE to change keywords
 * Aug  5 1996	Add HPUTNR8 to save specific number of decimal places
 * Oct 15 1996	Fix spelling
 * Nov  1 1996	Add DEG2STR to set specific number of decimal places
 * Nov  1 1996	Allow DEC2STR to handle upt to 6 decimal places
 *
 * Mar 20 1997	Fix format error in DEG2STR
 * Jul  7 1997	Fix 2 errors in HPUTCOM found by Allan Brighton
 * Jul 16 1997	Fix error in HPUTC found by Allan Brighton
 * Jul 17 1997	Fix error in HPUTC found by Allan Brighton
 */
