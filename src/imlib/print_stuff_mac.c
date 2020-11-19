/**
 ** functions for printing images on the laser_writer
 ** 	int		print_im(short **f, int N);	
 ** 	void	ps(char *string);
 ** 	void	PrPicDoc(PicHandle thePic, THPrint hPrint);
 **
 ** print_im() requires an 8-bit image (i.e. f-values of 0-256)
 **/
 
 
/*#define	MAC	*/		/* comment this line out for unix */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <signal.h>

#include "print_stuff.h"
#include "error.h"
#include "fits.h"

#define	PAGE_WIDTH			512
#define	PAGE_HEIGHT			700

int			psstringlen;
char		*psstring;
void		ps(char *string);
void		psDrawChar(char aChar);
void		print_caption(char *caption);
#ifdef MAC
#include <Quickdraw.h>
#include <Types.h>
#include <PrintTraps.h>
#include <Memory.h>
#include <OSUtils.h>
#define POSTSCRIPTBEGIN		190
#define POSTSCRIPTEND		191
#define POSTSCRIPTHANDLE	192
Handle		myHandle;
GrafPtr		MPWPort; 
void	PrPicDoc(PicHandle thePic, THPrint hPrint);
void	mac_open();
void	mac_close();
PicHandle 	pic;
GrafPort	plotPort;
#else
FILE		*lprpipe;
#endif

void		print_im(int N, int fmin, int fmax, char *caption)	
{
	char		pixval[256][2];
	char		hexchar[16] = {'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};
	int		i, j, fval;
	short		*f;
	extern int	ascii;
	
	for (i = 0; i < 16; i++) {
		for (j = 0; j < 16; j++) {
			pixval[255 - (16 * i + j)][0] = hexchar[i];
			pixval[255 - (16 * i + j)][1] = hexchar[j];
		}
	}

#ifndef MAC
/*	if ((lprpipe = fopen("temp.ps", "w")) == NULL) */
/*	if ((lprpipe = popen("lpr", "w")) == NULL)
		error_exit("can't open pipe to lpr\n"); */
#endif
	
	psstringlen = (2 * N + 2 > 1024 ? 2 * N + 2 : 1024);
	psstring = (char *) calloc(psstringlen, sizeof(char));

#ifdef	MAC	
	mac_open();
#else
	ps("\%!PS");
	sprintf(psstring, "0 %d translate", PAGE_HEIGHT);
	ps(psstring);
	ps("1 -1 scale");
#endif

	print_caption(caption);

	sprintf(psstring, "/picstr %d string def", N);
	ps(psstring);
	ps("0 0 translate");
	
	sprintf(psstring, "%d %d scale", PAGE_WIDTH, PAGE_WIDTH);
	ps(psstring);
	sprintf(psstring, 
		"/drawImage {%d %d 8 [%d 0 0 %d 0 0] {currentfile picstr readhexstring pop} image} def", 
		N, N, N, N);
	ps(psstring);
	ps("drawImage");
	f = (short *) calloc(N, sizeof(short));
	if (!f)
		error_exit("print_im: memory allocation failure\n");
	for (i = 0; i < N; i++) {
		if (ascii)
			read_ascii_line(f, N);
		else
			read_fits_line(f, N);
		for (j = 0; j < N; j++) {
			fval = 256 * (f[j] - fmin) / (float) (fmax - fmin);	
			fval = (fval >= 0 ? (fval < 256 ? fval : 255) : 0);
			psstring[2 * (N - j - 1)] = pixval[fval][0];
			psstring[2 * (N - j - 1) + 1] = pixval[fval][1];
		}
		psstring[2 * j] = '\0';
		ps(psstring);
	}
	free(f);

#ifdef	MAC	
	mac_close();
#else
	ps("showpage");
#endif
}

void	ps(char *string)
{
#ifdef	MAC
	long	len;

	HLock(myHandle);
	strcpy((char *) *myHandle, string);
	strcat((char *) *myHandle, "\n");				/* was slash-r */
	len = (long) strlen((char *) *myHandle);
	PicComment(192, len, myHandle);
	HUnlock(myHandle);
#else
/*	fprintf(lprpipe, "%s\n", string); */
	fprintf(stdout, "%s\n", string);
#endif
}


void	psDrawChar(char aChar)
{
	switch (aChar) {
		case '(':
			sprintf(psstring, "(\\50) hshow");
			break;
		case ')':
			sprintf(psstring, "(\\51) hshow");
			break;
		default:
			sprintf(psstring, "(%c) hshow", aChar);
			break;
	}
	ps(psstring);
}



void	print_caption(char *caption) {
	float	size = 10, space = 13, h = 13, v = 530;
	int		c, len;
	len = strlen(caption);
	
	ps("/hshow {1 -1 scale show 1 -1 scale} def"); 
	ps("/Times-Roman findfont   10.000 scalefont setfont");
	sprintf(psstring, "%f %f moveto", h, v);
	ps(psstring);
	for (c = 0; c < len; c++)
		switch (caption[c]) {
			case '\n':
				v += space;
				sprintf(psstring, "%f %f moveto", h, v);
				ps(psstring);
				break;
			default:
				psDrawChar(caption[c]);
		}
}

#ifdef	MAC
void	mac_open()
{
	Rect		pageRect;
	
	InitGraf((Ptr) &qd.thePort);				/* initialisation */
	GetPort(&MPWPort);
	OpenPort(&plotPort);
	SetRect(&pageRect, 0, 0, PAGE_WIDTH, PAGE_HEIGHT);
	ClipRect(&pageRect);
	pic = OpenPicture(&pageRect);
	PenSize(0, 0);
	MoveTo(0, 0);					/* seems to be necessary to print something */
	LineTo(0, 0);					/* to set the printer in a nice state		*/
	PenNormal();					/* god only knows why - just do it			*/
	PicComment(POSTSCRIPTBEGIN, 0, NULL);		/* from now on the laser writer */
	myHandle = NewHandle(psstringlen + 2);	
}


void	mac_close()
{
	THPrint		hPrint = 0;
	GrafPtr		savePort; 
	TPrStatus	prStatus;

	PicComment(POSTSCRIPTEND, 0, NULL);
	ClosePicture();
	SetPort(MPWPort);
    PrOpen();
	PrintDefault(hPrint = (TPrint **) NewHandle( sizeof( TPrint )));
	SetCursor( &qd.arrow );
	if (PrJobDialog(hPrint) != 0) {
		GetPort(&savePort);
		PrPicDoc (pic, hPrint);
		PrPicFile( hPrint, 0L, 0L, 0L, &prStatus );
		SetPort(savePort);
	}
	PrClose();
}


void	PrPicDoc(PicHandle thePic, THPrint hPrint)
{
	Rect 			printRect, r;
	TPPrPort		printPort;

	SetRect(&r, 0, 0, PAGE_WIDTH, PAGE_HEIGHT);
	printPort = PrOpenDoc( hPrint, 0L, 0L );
	SetPort((GrafPtr) printPort);
	printRect = (**hPrint).prInfo.rPage;
	PrOpenPage( printPort, 0L );
	DrawPicture(thePic, &r);
	PrClosePage( printPort );
	PrCloseDoc( printPort );
}
#endif



