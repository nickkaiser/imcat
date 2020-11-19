#include <stdio.h>
#define  CDL_LIBRARY_SOURCE
#include "cdl.h"


/*
 *  COMMUNICATIONS INTERFACE -- The communications interface handles all the
 *  low-level communications with the server.  It implements only the subset
 *  of the IIS protocol used by XImtool, SAOtng and SAOimage, not the entire
 *  IIS protocol.  It may be swapped out for another protocol in the future
 *  without affecting the tasks which use it as long as the basic steps re-
 *  quired for image display are the same.
 *  
 *            com_writeData  (fd, x, y, pix, nx, ny)
 *             com_readData  (fd, x, y, &pix, &nx, &ny)
 *           com_readCursor  (fd, sample, &x, &y, &key)
 *            com_setCursor  (fd, x, y, wcs)
 *          com_setFBConfig  (fd, configno)
 *             com_setFrame  (fd, frame)
 *             com_writeWCS  (fd, name, a, b, c, d, tx, ty, z1, z2, zt)
 *              com_readWCS  (fd, &name, &a, &b, &c, &d, &tx, &ty,
 *  				  &z1, &z2, &zt)
 *           com_eraseFrame  (fd)
 * 
 *	We do not actually display images here, all we do is send the individual
 *  set frame, write a data chunk, etc commands.  The caller is responsible for
 *  coordinating these calls into a legal sequence for image display.
 *	All routines return 0 if successfull and 1 otherwise.
 */

/* Function prototypes */
#ifdef __STDC__
#include <stddef.h>
#include <stdlib.h>
#endif


#define MEMORY          01              /* frame buffer i/o       */
#define LUT             02              /* lut i/o                */
#define FEEDBACK        05              /* used for frame clears  */
#define IMCURSOR        020             /* logical image cursor   */
#define WCS             021             /* used to set WCS        */
#define	SZ_IMCURVAL	160		/* cursor value str length*/
#define	SZ_WCSBUF	320		/* wcs buffer length      */

/* Command definitions */
#define PACKED          0040000
#define IMC_SAMPLE   	0040000
#define COMMAND         0100000
#define IIS_WRITE       0400000
#define IIS_READ        0100000

#define	OK		0
#define ERR		1

#undef	min
#define	min(a,b)	(a < b ? a : b)

/* IIS header packet structure.  DO NOT CHANGE. */
typedef struct  {
        short   tid;          		/* transfer id 	*/
        short   thingct;          	/* thing count 	*/
        short   subunit;          	/* subunit 	*/
        short   checksum;          	/* check sum 	*/
        short   x, y, z, t;		/* registers 	*/
} iis_hdr;   


static int 	frame = 1, fbconfig = 1;
static int	com_debug = 0;

#ifdef ANSI_FUNC

static int com_whdr(int fd, int tid, int subunit, int thingct, int x, int y, int z, int t);
static int com_write(int fd, char *buf, int nbytes);
static int com_read(int fd, char *buf, int maxbytes, int *nbytes);

#else

static int 	com_whdr(), com_write(), com_read();

#endif



/*  COM_WRITEDATA --  Write a block of data to the display.  This does not
 *  display the entire image, it just writes an array of pixels to the server.
 */

#ifdef ANSI_FUNC

int 
com_writeData (
    int fd,				/* connection file descriptor  	*/
    short x,
    short y,				/* corner of array		*/
    uchar *pix,				/* pixel array			*/
    int nbytes				/* number of bytes write		*/
)
#else

int
com_writeData (fd, x, y, pix, nbytes)
int	fd;				/* connection file descriptor  	*/
short	x, y;				/* corner of array		*/
uchar 	*pix;				/* pixel array			*/
int	nbytes;				/* number of bytes write		*/
#endif
{
    	/*  Send the IIS command for a data write. */
	if (com_whdr (fd, IIS_WRITE | PACKED, MEMORY, -nbytes, x, y,
	    1<<(frame-1), 0))
	        return (ERR);

        /* Send the pixels. */
        return (com_write (fd, (char *)pix, nbytes));
}


/*  COM_READDATA --  Read a block of data at a given location from the
 *  display.
 */

#ifdef ANSI_FUNC

int 
com_readData (
    int fdin,
    int fdout,			/* connection file descriptors 	*/
    short x,
    short y,				/* corner of readout		*/
    uchar *pix,				/* output pixel array		*/
    int *npix				/* number of bytes read		*/
)
#else

int
com_readData (fdin, fdout, x, y, pix, npix)
int	fdin, fdout;			/* connection file descriptors 	*/
short	x, y;				/* corner of readout		*/
uchar 	*pix;				/* output pixel array		*/
int	*npix;				/* number of bytes read		*/
#endif
{
	int nb = 0, n = *npix;

    	/*  Send the IIS command for a data read. */
	if (com_whdr (fdout, IIS_READ | PACKED, MEMORY, -(*npix), x, y,
	    1<<(frame-1), 0))
	        return (ERR);

        /* Get the pixels. */
	if (pix == NULL)
	    pix = (uchar *) malloc ((unsigned) *npix);

	while (nb < *npix) {
            if (com_read (fdin, (char *)pix+nb, n, &n))
		return (ERR);
	    nb += n;
	    n = *npix - nb;
	}
        return (OK);
}


/*  COM_READCURSOR --  Read the current cursor position.  If sample is set the
 *  value of the cursor is returned immediately, otherwise the read
 *  is blocked until the user hits a key.
 */

#ifdef ANSI_FUNC

int 
com_readCursor (
    int fdin,
    int fdout,			/* connection file descriptors 	*/
    int sample,				/* sample cursor or block	*/
    float *x,
    float *y,				/* output cursor coords		*/
    char *key				/* keystroke hit		*/
)
#else

int
com_readCursor (fdin, fdout, sample, x, y, key)
int	fdin, fdout;			/* connection file descriptors 	*/
int	sample;				/* sample cursor or block	*/
float	*x, *y;				/* output cursor coords		*/
char	*key;				/* keystroke hit		*/
#endif
{
	char	buf[SZ_IMCURVAL];
	int	n = SZ_IMCURVAL, wcs;

    	/* Send the IIS command for a cursor read. */
	if (com_whdr (fdout, (IIS_READ | (sample ? IMC_SAMPLE : 0)), IMCURSOR,
	    0, 0, 0, 0, 0))
		return (ERR);

    	/* Read back the ascii string and extract the cursor position.  */
	buf[0] = '\0';
    	if (com_read (fdin, buf, n, &n))
            return (ERR);

	*key = '\0';
    	return ((sscanf (buf, "%f %f %d %c", x, y, &wcs, key) != 4));
}


/*  COM_SETCURSOR -- Set the image cursor position.
 */

#ifdef ANSI_FUNC

int 
com_setCursor (
    int fd,				/* connection file descriptor  	*/
    int x,
    int y,				/* cursor coords		*/
    int wcs				/* cursor wcs			*/
)
#else

int
com_setCursor (fd, x, y, wcs)
int	fd;				/* connection file descriptor  	*/
int	x, y;				/* cursor coords		*/
int	wcs;				/* cursor wcs			*/
#endif
{
    	/* Pack up IIS command to set cursor */
	if (com_whdr (fd, IIS_WRITE, IMCURSOR, 0, x, y, wcs, 0))
	    return (ERR);
	return (OK);
}


/*  COM_SETFBCONFIG -- Set the frame buffer configuration number.  This is
 *  passed to the server in a WCS set command so we'll just save it for now.
 */

#ifdef ANSI_FUNC

int 
com_setFBConfig (
    int fd,				/* connection file descriptor  	*/
    int configno			/* fb configuration number	*/
)
#else

int
com_setFBConfig (fd, configno)
int	fd;				/* connection file descriptor  	*/
int	configno;			/* fb configuration number	*/
#endif
{
	fbconfig = configno;
	return (OK);
}


/*  COM_SETFRAME -- Set the current display frame.
 */

#ifdef ANSI_FUNC

int 
com_setFrame (
    int fd,				/* connection file descriptor  	*/
    int frame_num			/* frame number			*/
)
#else

int
com_setFrame (fd, frame_num)
int	fd;				/* connection file descriptor  	*/
int	frame_num;			/* frame number			*/
#endif
{
	short 	fr = 0;

	frame = frame_num;
	fr = 1 << (frame_num-1);

        /* Send the IIS command to select a frame.  */
	if (com_whdr (fd, IIS_WRITE, LUT | COMMAND, -1, 0, 0, 0, 0))
	    return (ERR);

        /* Send the frame info. */
        return (com_write (fd, (char *)&fr, sizeof(short)));
}


/*  COM_WRITEWCS --  Send the (linear) WCS to the server.
 */

#ifdef ANSI_FUNC

int 
com_writeWCS (
    int fd,				/* connection file descriptor  	*/
    char *name,				/* image name			*/
    float a,
    float b,
    float c,
    float d,			/* wcs coefficients		*/
    float tx,
    float ty,				/* translation values		*/
    float z1,
    float z2,				/* min/max			*/
    int zt				/* Z transform type		*/
)
#else

int
com_writeWCS (fd, name, a, b, c, d, tx, ty, z1, z2, zt)
int	fd;				/* connection file descriptor  	*/
char	*name;				/* image name			*/
float	a, b, c, d;			/* wcs coefficients		*/
float	tx, ty;				/* translation values		*/
float	z1, z2;				/* min/max			*/
int	zt;				/* Z transform type		*/
#endif
{
	char wcs_info[SZ_WCSBUF];
	int  nbytes;

    	/* Format the string to set world coordinate parameters.  */
    	(void) sprintf (wcs_info, "%s\n%g %g %g %g %g %g %g %g %d\n",
                name, a, b, c, d, tx, ty, z1, z2, zt);
	nbytes = strlen (wcs_info) + 1;

    	/* Send the IIS header to set a WCS. */
	if (com_whdr (fd, IIS_WRITE | PACKED, WCS, -nbytes, 0, 0, 1<<(frame-1),
	    fbconfig-1))
	        return (ERR);

    	/* Send the wcs info. */
    	return (com_write (fd, wcs_info, nbytes));
}


/*  COM_READWCS --  Read the (linear) WCS from the server.
 */

#ifdef ANSI_FUNC

int 
com_readWCS (
    int fdin,
    int fdout,			/* connection file descriptors 	*/
    char *name,				/* image name			*/
    float *a,
    float *b,
    float *c,
    float *d,			/* wcs coefficients		*/
    float *tx,
    float *ty,			/* translation values		*/
    float *z1,
    float *z2,			/* min/max			*/
    int *zt				/* Z transform type		*/
)
#else

int
com_readWCS (fdin, fdout, name, a, b, c, d, tx, ty, z1, z2, zt)
int	fdin, fdout;			/* connection file descriptors 	*/
char	*name;				/* image name			*/
float	*a, *b, *c, *d;			/* wcs coefficients		*/
float	*tx, *ty;			/* translation values		*/
float	*z1, *z2;			/* min/max			*/
int	*zt;				/* Z transform type		*/
#endif
{
	char wcs_info[SZ_WCSBUF];
	int	n = SZ_WCSBUF, tokens;

	/* Send the IIS command to read the WCS. */
	if (com_whdr (fdout, IIS_READ, WCS, 0, 0, 0, 0, 0))
	    return (ERR);

    	if (com_read (fdin, wcs_info, n, &n))
            return (ERR);

    	tokens = sscanf (wcs_info, "%[^\n]\n%g%g%g%g%g%g%g%g%d",
	    name, a, b, c, d, tx, ty, z1, z2, zt);
    	if (tokens == EOF) {
            /* Had no WCS info, can't even get name.  no error.  */
            name[0] = '\0';
            *a = *b = *c = *d = *tx = *ty = *z1 = *z2 = 0.0;
            *zt = 0;

    	} else if (tokens < 10)
            /* partial read, something must be wrong.  */
            return (ERR);

	return (OK);
}


/*  COM_ERASEFRAME -- Clear the current frame.
 */

#ifdef ANSI_FUNC

int 
com_eraseFrame (
    int fd				/* connection file descriptor 	*/
)
#else

int
com_eraseFrame (fd)
int	fd;				/* connection file descriptor 	*/
#endif
{
    	/* Send IIS command to erase. */
	if (com_whdr (fd, IIS_WRITE+(fbconfig-1), FEEDBACK, 0, 0, 0,
	    1<<(frame-1), 0))
	    return (ERR);
	return (OK);
}



/*------------------
  PRIVATE PROCEDURES
  ------------------*/


/* COM_WHDR -- Load and send the IIS header packet given the structure elements.
 */

#ifdef ANSI_FUNC

static int 
com_whdr (
    int fd,				/* connection file descriptor 	*/
    int tid,				/* thing ID			*/
    int subunit,			/* subunit			*/
    int thingct,			/* thing count			*/
    int x,
    int y,
    int z,
    int t			/* registers			*/
)
#else

static int
com_whdr (fd, tid, subunit, thingct, x, y, z, t)
int	fd;				/* connection file descriptor 	*/
int	tid;				/* thing ID			*/
int	subunit;			/* subunit			*/
int	thingct;			/* thing count			*/
int	x, y, z, t;			/* registers			*/
#endif
{
	iis_hdr	iis;
    	int sum = 0;

	/* Load the structure. */
    	iis.tid      = (short) tid;
    	iis.subunit  = (short) subunit;
    	iis.thingct  = (short) thingct;
    	iis.checksum = (short) 0;
    	iis.x 	     = (short) x;
    	iis.y 	     = (short) y;
    	iis.z 	     = (short) z;
    	iis.t 	     = (short) t;

	/* Compute the checksum. */
        sum = iis.tid + iis.subunit + iis.thingct + iis.checksum +
    	      iis.x + iis.y + iis.z + iis.t;
    	iis.checksum = 0177777 - sum;

	if (com_debug) {
            printf (
               "subunit=%06o tid=%06o nbytes=%7d x=%05o y=%05o z=%05o t=%05o\n",
               iis.subunit & 077,
               iis.tid,
	       (!(iis.tid & PACKED) ? (-iis.thingct * 2) : (-iis.thingct)),
               iis.x & 0177777,
               iis.y & 0177777,
               iis.z & 0177777,
               iis.t & 0177777);
            fflush (stdout);
	}

    	/* Send the header and return the success code.  */
    	return (com_write(fd, (char *)&iis, sizeof(iis)));
}


/* COM_WRITE -- Asynchronous write of data to the server.  Write exactly
 * nbytes bytes from the buffer to the server.
 */

#ifdef ANSI_FUNC

static int 
com_write (
    int fd,				/* connection file descriptor 	*/
    char *buf,				/* buffer to write		*/
    int nbytes				/* number of bytes to write	*/
)
#else

static int
com_write (fd, buf, nbytes)
int	fd;				/* connection file descriptor 	*/
char 	*buf;				/* buffer to write		*/
int 	nbytes;				/* number of bytes to write	*/
#endif
{
    	int n = 0, total = 0, maxbytes = nbytes;
	char *ip = (char *)buf;

	for (total=0; total < nbytes; total += n, ip += n) {
	    n = nbytes - total;
	    if (maxbytes)
		n = min (maxbytes, n);
            if ((n = write (fd, ip, n)) <= 0)
                return (ERR);
	}
    	return (OK);
}


/* COM_READ -- Read data from the server.  Try to read at most maxbytes bytes
 * from the server into the buffer, return the number of bytes actually read.
 */

#ifdef ANSI_FUNC

static int 
com_read (
    int fd,				/* connection file descriptor 	*/
    char *buf,				/* buffer to write		*/
    int maxbytes,			/* max number of bytes to read	*/
    int *nbytes			/* number of bytes actually read*/
)
#else

static int
com_read (fd, buf, maxbytes, nbytes)
int	fd;				/* connection file descriptor 	*/
char 	*buf;				/* buffer to write		*/
int 	maxbytes;			/* max number of bytes to read	*/
int 	*nbytes;			/* number of bytes actually read*/
#endif
{
    	if ((*nbytes = read (fd, buf, maxbytes)) <= 0) {
            buf[0] = '\0';
            *nbytes = 0;
            return (ERR);
    	}
    	return (OK);
}


/* COM_SETDEBUG -- Set the state of the debug flag.
 */

#ifdef ANSI_FUNC

int 
com_setDebug (int state)
#else

int
com_setDebug (state)
int     state;
#endif
{
        com_debug = state;
        return (OK);
}
