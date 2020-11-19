/*
 * mosaicmap2.c
 *
 * maps a source fits file to a target fits file using transformations
 * computed by mosaicfit2
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../utils/arrays.h"
#include "../../imlib/fits.h"
#include "../../imlib/map.h"
#include "../../utils/modefunc.h"

#define usage  "\n\n\
NAME\n\
	mosaicmap2 --- maps fits file using transformations computed by mosaicfit2\n\
\n\
SYNOPSIS\n\
	mosaicmap2 alpha Mx My x0 y0 mx my N1 N2 parfile srcfits dstfits mapmode\n\
\n\
DESCRIPTION\n\
	'mosaicmap2' maps a source fits file to a target file using transformation computed by\n\
	'mosaicfit2'. Alpha is the assumed distortion model parameter,\n\
	Mx, My define the margins by which the target image overlaps the source image,\n\
	x0, y0 is the location of pottom left pixel in nominal coords,\n\
	mx, my are the margin parameters as defined in 'nominal.db',\n\
	N1, N2 are such that the target image size is N1 + 2 * Mx, N2 + 2 * My\n\
	(these would normally be the dimensions of the source image)\n\
	and 'parfile' contains the transformation parameters.\n\
	Srcfits, dstfits are the source and target image files.\n\
	If mapmode is 0, 1, or 2, then mosaicmap2 will actually map the source image\n\
	using nearest pixel, linear interpolation and triangular\n\
	tesselation mapping modes respctively.\n\
	If mapmode is a negative integer -n, then it generates a 2^(n-1) times\n\
	scrunched deflection image, which can then be fed to 'mapbynumericdef'\n\
	to actually do the mapping. In this mode srcfits is ignored.\n\
	Thus, for n = -1 we generate a full size deflection image, n = -4\n\
	gives factor 8 linear size compression etc. In this mode, the\n\
	source image name is ignorred.\n\
	'mosaicmap2' is meant to be invoked by ''mosaicfit2'.\n\
	Give parfile = 'NULL' to map the reference image.\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n"

/* global parameters - needed by deflection() */
double		alpha, Mx, My, xo, yo, mx, my;
int             nmodes, *l, *m;
double           **a;

void	deflect(float ri, float rj, float *di, float *dj);

main(int argc, char *argv[]) {
	int	N1, N2, NX, NY, NXsrc, NYsrc, comc, mapmode, makedefim, i, ix, iy, scalefac, scrunchfac; 
	char 	*comv[MAX_COMMENTS], *parfilename;
	float   **fsource, **ftarget, **defx, **defy;
	FILE 	*ipf, *opf, *parf;
	int	ndim;
	char	*vardef[MODEFUNC_MAX_VARS], *xname;
	int	nvar;

	/* defaults */
	makedefim = 0;

	if (argc != 14) {
		fprintf(stderr, usage);
		exit(-1);
	}
	sscanf(argv[1], "%lf", &alpha);
	sscanf(argv[2], "%lf", &Mx);
	sscanf(argv[3], "%lf", &My);
	sscanf(argv[4], "%lf", &xo);
	sscanf(argv[5], "%lf", &yo);
	sscanf(argv[6], "%lf", &mx);
	sscanf(argv[7], "%lf", &my);
	sscanf(argv[8], "%d", &N1);
	sscanf(argv[9], "%d", &N2);
	parfilename = argv[10];
	sscanf(argv[13], "%d", &mapmode);
	if (mapmode < 0) {
		scrunchfac = - 1 - mapmode;
		makedefim = 1;
		scalefac = 1;
		for (i = 0; i < scrunchfac; i++) {
			scalefac *= 2;
		}
	}

	/* read the source image if necessary */
	if (!makedefim) {
		if (!(ipf = fopen(argv[11], "r"))) {
			error_exit("mosaicmap2: unable to open source image file\n");
		}
		set_fits_ipf(ipf);
		fread_fits(&fsource, &NXsrc, &NYsrc, &comc, comv);
		fclose(ipf);
	}

	/* check that we can actually open the target image */
	if (!(opf = fopen(argv[12], "w"))) {
		error_exit("mosaicmap2: unable to open target image file\n");
	}

	/* get the mapping parameters */
	if (strcmp(parfilename, "NULL")) {
		get2Dpolymodel(parfilename, &l, &m, &ndim, &a, &nmodes, &nvar, vardef, &xname);
		if (ndim != 2) {
			error_exit("warpimage: I need a 2D parameter file!\n");
		}
	} else {
		nmodes = 0;
	}

	/* compute output image size */
	NX = N1 + 2 * Mx;
	NY = N2 + 2 * My;

	/* now M becomes margin relative to nominal tile size */
	Mx -= mx;
	My -= my;

	if (makedefim) {
		NX = (int) floor(NX / scalefac);
		NY = (int) floor(NY / scalefac);
		allocFloatArray(&ftarget, NX, 2 * NY);
		defx = ftarget;
		defy = ftarget + NY;
		for (iy = 0; iy < NY; iy++) {
			for (ix = 0; ix < NX; ix++) {
				deflect((float) (iy * scalefac), (float) (ix * scalefac), &(defy[iy][ix]), &(defx[iy][ix]));
			} 
		}
	} else {
		allocFloatArray(&ftarget, NX, NY);
		switch (mapmode) {
			case ULTRAFAST_MAP_MODE:
				ultrafastmap(ftarget, NX, NY, fsource, NXsrc, NYsrc, deflect);
				break;
			case FAST_MAP_MODE:
				fastmap(ftarget, NX, NY, fsource, NXsrc, NYsrc, deflect);
				break;
			case TRIANGLE_MAP_MODE:
				map(ftarget, NX, NY, fsource, NXsrc, NYsrc, deflect);
				break;
			default:
				fprintf(stderr, "mosaicmap: bad mapping mode\n");
				exit(-1);
				break;
		}
	}

	/* write the target image */
	set_fits_opf(opf);
        add_comment(argc, argv, &comc, comv);
	if (makedefim) {
		set_output_pixtype(FLOAT_PIXTYPE);
        	fwrite_fits(ftarget, NX, 2 * NY, comc, comv);
	} else {
        	fwrite_fits(ftarget, NX, NY, comc, comv);
	}

        /* all done */
        exit(0);
	
}


void	deflect(float ri, float rj, float *di, float *dj)
{
        int     i, mode;
        double  fmode, ff[2], r[2], d[2], rp[2], rr;

        r[0] = (double) (rj + xo - (mx + Mx));
        r[1] = (double) (ri + yo - (my + My));
        for (i = 0; i < 2; i++) {
                ff[i] = 0.0;
        }
        for (mode = 0; mode < nmodes; mode++) {
                fmode = f(l[mode], m[mode], r);
                for (i = 0; i < 2; i++) {
                        ff[i] += a[i][mode] * fmode;
                }
        }
	rp[0] = r[0] + ff[0];
	rp[1] = r[1] + ff[1];
	rr = rp[0] * rp[0] + rp[1] * rp[1];
	d[0] = rp[0] * alpha * rr;
	d[1] = rp[1] * alpha * rr;
        *dj = (float) (ff[0] + d[0] - Mx - mx);
        *di = (float) (ff[1] + d[1] - My - my);
}

