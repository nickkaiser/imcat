/*
 * iis.c
 *
 * create iis header followed by an image line
 *
 */

#define usage "\n\n\
NAME\n\
	iis - pipe fits image data into saopipe or ximtool in IIS format\n\
\n\
SYNOPSIS\n\
	iis [-u] [-f fmin fmax] [-F fitsfile] [-p pscom] [-P pscom]\n\
		[-c catfile [-s size] [-t type] [-T tname] [-e ename]\n\
		[-r rname] [-n nname] [-l]]\n\
\n\
DESCRIPTION\n\
	\"iis\" reads a fits image from stdin, linearly scales the\n\
	pixel values to range 0-199, prepends an iis header\n\
	and writes the output to a FIFO pipe /dev/imt1o so that saoimage\n\
	or ximtool can display it. This saves all the fiddle faddle\n\
	with the viewers GUI and is useful to get visual feedback\n\
	from scripts which process a series of images.\n\
\n\
	The first version of this program was devised from a fax of\n\
	a photocopy of the arcane IIS display device protocol from\n\
	George Miyashiro.  Karl Glazebrook then refined this, improving\n\
	efficiency - thanks Karl - but I have now ditched his library\n\
	in favour of the IRAF/NOAO cdl library due to Michael\n\
	Fitzpatrick, which does the same thing, but which is\n\
	better documented and has some useful extensions.\n\
\n\
	Bad pixels (flagged by the MAGIC value) are highlighted in green.\n\
\n\
	Options:\n\
\n\
	-u	print this mesage\n\
	-f	set the limits for pixel values\n\
	-F	read the image from 'fitsfile'\n\
	-p	send postscript output to command 'pscom'\n\
	-P	send annotated postscript output to command 'pscom'\n\
	-c 	read in the objects from 'catfile' and display them\n\
		as overlay markers. If catfile = '-' the catalogue\n\
		will be read from stdin (but be sure then to use -F\n\
		option also). With -c you can use the following options\n\
		to control the appearance of markers:\n\
	-s 	set default size of markers (5)\n\
	-t	set default type (8). Allowed types are:
			POINT		2\n\
			BOX		4\n\
			PLUS		8\n\
			CROSS		16\n\
			DIAMOND		32\n\
			CIRCLE		64\n\
			STAR		128\n\
			HLINE		256\n\
			VLINE		512\n\
		you can add these to make composite markers and/or\n\
		add 1 to fill them.\n\
	-e 	option to draw ellipses with polarisation ename\n\
	-r 	set size to input variable rname\n\
	-T	set type from inut variable tname\n\
	-n	label marker with input variable 'nname'\n\
	-C	set color from input variable 'cname'. Allowed colors are\n\
			BLACK	1	CORAL	9\n\
			WHITE	2	MAROON	10\n\
			RED	3	ORANGE	11\n\
			GREEN	4	KHAKI	12\n\
			BLUE	5	ORCHID	13\n\
			YELLOW	6	TQUOISE	14\n\
			CYAN	7	VIOLET	15\n\
			MAGENTA	8	WHEAT	16\n\
	-l	enter edit loop after displaying cat in which we do the following\n\
		in response to key strokes:\n\
			Key		Action\n\
			---		---\n\
			d		delete marker nearest cursor\n\
			q		quit edit loop\n\
			[0-9]		set overlay color\n\
			l		start or end a line\n\
			b		start or end a box\n\
			r		redraw overlay\n\
			f		toggle fill state for boxes\n\
\n\
AUTHOR\n\
	Nick Kaiser -- kaiser@hawaii.edu\n\
\n\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <values.h>
#include <unistd.h>
#include "../../imlib/fits.h"
#include "../../utils/error.h"
#include "../../cdllib/cdl.h"
#include "../../cdllib/cdlP.h"
#include "../../catlib/cat.h"

#define MAGIC		FLOAT_MAGIC
#define	PIX_MAX		200
#define	DIM_MAX		4096
#define	COLORBASE	201
#define	NKBOX		0
#define	NKLINE		1

main(int argc, char *argv[]) {
	float	**f, fmin, fmax, fscaled, tx, ty;
	int	ix, iy, N1, N2, arg = 1, frame, autoscale, markcat;
	fitsheader	*fits;
	CDLPtr 	cdl;
	uchar	*pix;
	int	fb, w, h, nf, lx, ly, n1, n2;
	char	key[1], *catfilename, *imfilename, *xname, defaultxname[2] = "x", *pscom;
	float	xc, yc, xc1, yc1, xc2, yc2, phi, emod;
	int	nobj, a, b;
	extern	MarkerPtr	DLHead[], DLTail[];
	MarkerPtr		themarker;
	FILE	*catfile, *imfile;
	cathead	*ipcat;
	object	*ipobj, *ipobjbase;
	item	*xitem, *ritem, *eitem, *titem, *coloritem, *numberitem;
	int	xindex, rindex, eindex, tindex, colorindex, numberindex;
	double	*x, *r, *e, *t, *color, *number;
	char	*rname, *ename, *tname, *colorname, *numbername;
	int	variablesize, ellipses, defaultsize, defaulttype, variabletype, variablecolor;
	int	msize, mtype, defaultcolor, mcolor, numberlabels, mnumber, doloop, annotateps;
	int	fillboxes, catheadwritten;

	/* defaults */
	autoscale = 1;
	fmin = 0.0;
	fmax = 0.0;
	frame = 1;
	markcat = 0;
	xname = defaultxname;
	variablesize = 0;
	variabletype = 0;
	variablecolor = 0;
	ellipses = 0;
	defaultsize = 5;
	defaulttype = M_PLUS;
	defaultcolor = C_YELLOW;
	imfile = stdin;
	numberlabels = 0;
	doloop = 0;
	pscom = NULL;
	annotateps = 0;
	fillboxes = 0;
	catheadwritten = 0;

	/* parse args */
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			error_exit(usage);
		}
		switch (argv[arg++][1]) {
			case 'f':
				sscanf(argv[arg++], "%f", &fmin);
				sscanf(argv[arg++], "%f", &fmax);
				autoscale = 0;
				break;
			case 'F':
				imfilename = argv[arg++];
				imfile = fopen(imfilename, "r");
				if (!imfile) {
					error_exit("iis: failed to open image file\n");
				}
				break;
			case 'P':
				annotateps = 1;
			case 'p':
				pscom = argv[arg++];
				break;
			case 'c':
				markcat = 1;
				if (!strcmp(argv[arg], "-")) {
					catfile = stdin;
				} else {
					catfilename = argv[arg];
					catfile = fopen(catfilename, "r");
					if (!catfile) {
						error_exit("iis: failed to open cat file\n");
					}
				}
				arg++;
				break;	
			case 's':
				sscanf(argv[arg++], "%d", &defaultsize);
				break;
			case 't':
				sscanf(argv[arg++], "%d", &defaulttype);
				break;
			case 'T':
				variabletype = 1;
				tname = argv[arg++];
				break;
			case 'e':
				ellipses = 1;
				ename = argv[arg++];
				break;
			case 'r':
				variablesize = 1;
				rname = argv[arg++];
				break;
			case 'n':
				numberlabels = 1;
				numbername = argv[arg++];
				break;
			case 'C':
				variablecolor = 1;
				colorname = argv[arg++];
				break;
			case 'l':
				doloop = 1;
				break;
			case 'u':
			default:
				error_exit(usage);
				break;
		}
	}

	if ((imfile == stdin) && (catfile == stdin)) {
		error_exit("iis: can't read two standard inputs\n");
	}

	cdl = cdl_open("fifo:/dev/imt1i:/dev/imt1o");
	if (cdl == NULL) {
		error_exit("iis: failed to open fifo pipe to server\n");
	}

	cdl_displayFITS(cdl, "tmp.fits", 1, FB_AUTO, 0);
	cdl_close(cdl);
	exit(0);

	/* read image and reduce dimensions to 4096 if necessary */
	read2Dfloatimage(&f, &N1, &N2, &fits, imfile);
	N1 = (N1 > DIM_MAX ? DIM_MAX : N1);
	N2 = (N2 > DIM_MAX ? DIM_MAX : N2);
	
	/* do autoscaling if required */
	if (autoscale) {
		fmin = 1.e20;
		fmax = -1.e20;
		for (iy = 0; iy < N2; iy++) {
			for (ix = 0; ix < N1; ix++) {
				if (f[iy][ix] != (float) MAGIC) {
					if (f[iy][ix] > fmax)
						fmax = f[iy][ix];
					if (f[iy][ix] < fmin)
						fmin = f[iy][ix];
				}
			}
		}
	}
	if (fmin == fmax) {
		fmax += 1.0;
	}

	/* compute 8 bit image */
	pix = (uchar *) calloc(N1 * N2, sizeof(uchar));
	for (iy = 0; iy < N2; iy++) {
		for (ix = 0; ix < N1; ix++) {
			if (f[iy][ix] == MAGIC) {
				pix[N1 * iy + ix] = C_GREEN;
			} else {
				fscaled = PIX_MAX * (f[iy][ix] - fmin) / (fmax - fmin);
				fscaled = (fscaled < 0 ? 0 : (fscaled > PIX_MAX ? PIX_MAX : fscaled));
				pix[N1 * iy + ix] = (uchar) fscaled;
			}
		}
	}
	
	/* choose the frame buffer */
	n1 = (N1 >= 512 ? N1 : 512);
	n2 = (N2 >= 512 ? N2 : 512);
	cdl_selectFB(cdl, n1, n2, &fb, &w, &h, &nf, 1);
	fprintf(stderr, "# iis: frame buffer: fb=%d w=%d h=%d nf=%d\n", fb, w, h, nf);

	/* set the WCS */
	tx = (N1/2) - (w/2) + 0.5;
	/* slight idiosyncracy in the y WCS scaling */
	if (h == N2) {
		ty = (N2/2) + (h/2) - 0.5;
	} else {
		ty = (N2/2) + (h/2) + 0.5;
	}
/*	cdl_setWCS(cdl, "imcat", "iis", 1.0, 0.0, 0.0, -1.0, tx, ty, fmin, fmax, 1);*/

	/* set and clear the frame */
	cdl_setFrame(cdl, frame);
	cdl_clearFrame(cdl);

	/* compute offsets and display the image */
	lx = (w / 2) - (N1 / 2);
	ly = h - ((h/2) + (N2/2));
	cdl_writeSubRaster(cdl, lx, ly, N1, N2, pix);

	cdl_close(cdl);
	if (markcat) {
		/* read cat into linked list */
		setcatipf(catfile);
		ipcat = readcathead();
		xitem = getobjectitem(xname, ipcat);
		if ((xitem->itype != NUM_TYPE) || (xitem->ndim != 1) || ((xitem->dim)[0] != 2)) {
			error_exit("iis: position variable must be numerical 2-vector\n");
		}	
		ipobj = newobject(ipcat);
		xindex = getobjectitemindex(xname, ipobj);
		if (variabletype) {
			titem = getobjectitem(tname, ipcat);
			tindex = getobjectitemindex(tname, ipobj);
		}
		if (variablesize) {
			ritem = getobjectitem(rname, ipcat);
			rindex = getobjectitemindex(rname, ipobj);
		}
		if (variablecolor) {
			coloritem = getobjectitem(colorname, ipcat);
			colorindex = getobjectitemindex(colorname, ipobj);
		}
		if (numberlabels) {
			numberitem = getobjectitem(numbername, ipcat);
			numberindex = getobjectitemindex(numbername, ipobj);
		}
		if (ellipses) {
			eitem = getobjectitem(ename, ipcat);
			eindex = getobjectitemindex(ename, ipobj);
		}
		while (1) {
			ipobj = newobject(ipcat);
			allocobjectcontents(ipobj);
			if (!readobject(ipobj)) {
				break;
			}
			ipobj->next = ipobjbase;
			ipobjbase = ipobj;
			x = (double *) ((ipobj->addrlist)[xindex]);
			ix = (int) floor(x[0] + lx + 1);
			iy = (int) floor(x[1] + ly + 1);
			if (variablesize) {
				r = (double *) ((ipobj->addrlist)[rindex]);
				msize = (int) floor(0.5 + 2 * r[0]);
			} else {
				msize = 2 * defaultsize;
			}
			if (variabletype) {
				t = (double *) ((ipobj->addrlist)[tindex]);
				mtype = (int) floor(0.5 + t[0]);
			} else {
				mtype = defaulttype;
			}
			if (variablecolor) {
				color = (double *) ((ipobj->addrlist)[colorindex]);
				mcolor = (int) floor(COLORBASE + 0.5 + color[0]);
			} else {
				mcolor = defaultcolor;
			}
			if (numberlabels) {
				number = (double *) ((ipobj->addrlist)[numberindex]);
				mnumber = (int) floor(0.5 + number[0]);				
			} else {
				mnumber = 0;
			}
			if (ellipses) {
				e = (double *) ((ipobj->addrlist)[eindex]);
				phi = 90.0 * atan2(e[1], e[0]) / M_PI;
				emod = sqrt(e[0] * e[0] + e[1] * e[1]);
				if (variablesize) {
					a = (int) floor(0.5 + r[0] * (1 + emod));
					b = (int) floor(r[0] / (1 + emod));
				} else {
					a = (int) floor(defaultsize * (1 + emod));
					b = (int) floor(defaultsize / (1 + emod));
				}
				cdl_markEllipse(cdl, ix, iy, a, b, phi, 0, mcolor);
			} else {
				cdl_markPoint(cdl, ix, iy, mnumber, msize, mtype, mcolor);
			}
			DLTail[frame-1]->lcobject = ipobj;
		}
		if (doloop) {
			fprintf(stderr, "# iis: entering edit loop: active keys d,q,l,b,r,f,e,[0-9] \n");
		}
		while (doloop) {
			cdl_readCursor(cdl, 0, &xc1, &yc1, key);
			fprintf(stderr, "# iis: cursor event: %10.2f %10.2f %s\n", xc1, yc1, key);
			switch (key[0]) {
				case 'd':
					fprintf(stderr, "# iis: deleting object\n");
					cdl_deleteMark(cdl, (int) xc1 + lx, (int) yc1 + ly);
					break;
				case 'q':
					fprintf(stderr, "# iis: exiting edit loop\n");
					doloop = 0;
					break;
				case 'l':
					fprintf(stderr, "# iis: starting line:  hit 'l' to finish, any other key to exit\n");
					cdl_readCursor(cdl, 0, &xc2, &yc2, key);
					if (key[0] == 'l') {
						cdl_markLine(cdl, (int) xc1 + lx, (int) yc1 + ly, (int) xc2 + lx, (int) yc2 + ly, mcolor);
						DLTail[frame-1]->lcobject = NULL;
					}
					break;
				case 'b':
					fprintf(stderr, "# iis: starting box:  hit 'b' to finish, any other key to exit\n");
					cdl_readCursor(cdl, 0, &xc2, &yc2, key);
					if (key[0] == 'b') {
						cdl_markBox(cdl, (int) xc1 + lx, (int) yc1 + ly, (int) xc2 + lx, (int) yc2 + ly, fillboxes, mcolor);
						DLTail[frame-1]->lcobject = NULL;
					}
					break;
				case 'r':
					fprintf(stderr, "# iis: redrawing overlay\n");
					cdl_redrawOverlay(cdl);
					break;
				case 'f':
					fillboxes = (fillboxes ? 0 : 1);
					fprintf(stderr, "# iis: toggling fillboxes to %d\n", fillboxes);
					break;				
				case 'e':
					themarker = cdl_findNearest(DLHead[frame-1], (int) xc1 + lx, (int) yc1 + ly);
					if (!catheadwritten) {
						fprintf(stderr, "# iis: writing cathead to stdout\n");
						addargscomment(argc, argv, ipcat);
						setcatopfiletype(TEXT_FILE_TYPE);
						writecathead(ipcat);
						catheadwritten = 1;
					}
					if (themarker->lcobject != NULL) {
						fprintf(stderr, "# iis: writing object to stdout\n");
						writeobject(themarker->lcobject);
					}
					break;
				default:
					if (1 == sscanf(key, "%d", &mcolor)) {
						fprintf(stderr, "# iis: changing overlay color to %d\n", mcolor);
						mcolor += COLORBASE;
					} else {
						fprintf(stderr, "# iis: key %s not recognized\n", key);
					}
			}
		}
	}
	
	/* generate postsript */
	if (pscom) {
		fprintf(stderr, "# iis: printing postscript to command '%s'\n", pscom);
		cdl_readSubRaster(cdl, lx, ly, N1, N2, &pix);
		cdl_printPix(cdl, pscom, pix, N1, N2, annotateps);
	}
	cdl_close(cdl);
	exit(0);
}
