/* 
 * lmodel.c - functions to create and manipulate linear fit models ('lmodels')
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "catlib/cat.h"
#include "utils/error.h"
#include "lmodel.h"



/* write out a polymodel structure as lc-format cat */
int	writelmodel(lmodel *themodel, FILE *opstream)
{
	item 	*theitem, *pitem, *litem, *nitem, *mitem, *kitem, *iitem, *covaritem;
	int	dim, d, m;

	/* create the model type header item */
	theitem = newitem("model_type", TEXT_TYPE, 1, 1);
	allocitemcontents(theitem, &(theitem->addr), 0);
	((char **)(theitem->addr))[0] = (char *) calloc(128, sizeof(char));
	switch (themodel->modeltype) {
		case POLYNOMIAL_LMODEL:
			strcpy(((char **)(theitem->addr))[0], "polynomial");
			break;
		case ZERNIKE_LMODEL:
			strcpy(((char **)(theitem->addr))[0], "zernike");
			break;
		case FOURIER_LMODEL:
			strcpy(((char **)(theitem->addr))[0], "fourier");
			break;
		default:
			error_exit("writelmodel: invalid case\n");
	}
	installitem(theitem, &((themodel->cat)->headeritembase));

	/* and the other required header items */
	install1textheaderitem("aname", themodel->aname, themodel->cat);
	install1textheaderitem("xname", themodel->xname, themodel->cat);
	install1numericheaderitem("xdim", (double) (themodel->xdim), themodel->cat);

	/* optional header items */
	install1numericheaderitem("nmodes", (double) (themodel->nmodes), themodel->cat);
	if (themodel->hasxorigin) {
		theitem = newitem("xorigin", NUM_TYPE, 1, themodel->xdim);
		allocitemcontents(theitem, &(theitem->addr), 0);
		theitem->addr = (void *) (themodel->xorigin);
		installitem(theitem, &((themodel->cat)->headeritembase));
	}
	install1numericheaderitem("hascovar", themodel->hascovar, themodel->cat);

	switch (themodel->modeltype) {
		case POLYNOMIAL_LMODEL:
			install1numericheaderitem("lmin", (double) (themodel->lmin), themodel->cat);
			install1numericheaderitem("lmax", (double) (themodel->lmax), themodel->cat);
			break;
		case ZERNIKE_LMODEL:
			install1numericheaderitem("nmin", (double) (themodel->nmin), themodel->cat);
			install1numericheaderitem("nmax", (double) (themodel->nmax), themodel->cat);
			break;
		case FOURIER_LMODEL:
			install1numericheaderitem("kmin", (double) (themodel->kmin), themodel->cat);
			install1numericheaderitem("kmax", (double) (themodel->kmax), themodel->cat);
			install1numericheaderitem("lbox", (double) (themodel->lbox), themodel->cat);
			install1numericheaderitem("boxlimits", (double) (themodel->boxlimits), themodel->cat);
			break;
		default:
			error_exit("writelmodel: invalid case\n");
	}

	/* add the object items */
	switch (themodel->modeltype) {
		case POLYNOMIAL_LMODEL:
			pitem = newitem("p", NUM_TYPE, 1, themodel->xdim);
			addobjectitem(pitem, themodel->cat);
			litem = newitem("l", NUM_TYPE, 1, themodel->xdim);
			addobjectitem(litem, themodel->cat);
			break;
		case ZERNIKE_LMODEL:
			nitem = newitem("n", NUM_TYPE, 1, 1);
			addobjectitem(nitem, themodel->cat);
			mitem = newitem("m", NUM_TYPE, 1, 1);
			addobjectitem(mitem, themodel->cat);
			break;
		case FOURIER_LMODEL:
			kitem = newitem("k", NUM_TYPE, 1, themodel->xdim);
			addobjectitem(kitem, themodel->cat);
			iitem = newitem("i", NUM_TYPE, 1, 1);
			addobjectitem(iitem, themodel->cat);
			break;
		default:
			error_exit("writelmodel: invalid case\n");
	}
	addobjectitem(themodel->aitem, themodel->cat);
	if (themodel->hascovar) {
		covaritem = newitem("covar", NUM_TYPE, 1, themodel->nmodes);
		addobjectitem(covaritem, themodel->cat);
	}
	
	setcatopf(opstream);

	writecathead(themodel->cat);

	for (m = 0; m < themodel->nmodes; m++) {
		fprintf(opstream, " ");
		switch (themodel->modeltype) {
			case POLYNOMIAL_LMODEL:
				pitem->addr = (void *) ((themodel->p)[m]);
				writeitem(pitem, pitem->addr, 0);
				litem->addr = (void *) ((themodel->l)[m]);
				writeitem(litem, litem->addr, 0);
				break;
			case ZERNIKE_LMODEL:
				nitem->addr = (void *) &((themodel->n)[m]);
				writeitem(nitem, nitem->addr, 0);
				mitem->addr = (void *) &((themodel->m)[m]);
				writeitem(mitem, mitem->addr, 0);
				break;
			case FOURIER_LMODEL:
				kitem->addr = (void *) ((themodel->k)[m]);
				writeitem(kitem, kitem->addr, 0);
				iitem->addr = & ((themodel->i)[m]);
				writeitem(iitem, iitem->addr, 0);
				break;
			default:
				error_exit("writelmodel: invalid case\n");
		}
		writeitem(themodel->aitem, (themodel->a)[m], 0);
		if (themodel->hascovar) {
			writeitem(covaritem, (themodel->covar)[m], 0);
		}
		fprintf(opstream, "\n");
	}
			
	return(1);
}


lmodel	*readlmodel(FILE *ipstream)
{
	lmodel	*themodel;
	cathead	*ipcat;
	object	*ipobj, *ipobjbase;
	char	**strptr, *xname, *modeltypestring, *aname;
	double	*dblptr;
	item	*aitem, *xitem, *covaritem, *xoriginitem;
	int	ll, mm, nn, m, d, xdim, nmodes, modeltype;
	int	pindex, lindex, mindex, nindex, aindex, kindex, iindex, covarindex;

	setcatipf(ipstream);
	ipcat = readcathead();
	ipobj = newobject(ipcat);
	allocobjectcontents(ipobj);

	/* get the modeltype */
	strptr = (char **) getheaderitemaddress("model_type", ipcat);
	if (!strptr) error_exit("readlmodel : missing header item 'model_type'\n");
	modeltypestring = strptr[0];
	if (!strcmp(modeltypestring, "polynomial") || !strcmp(strptr[0], "2Dpoly")) {
		modeltype = POLYNOMIAL_LMODEL;
	} else if (!strcmp(modeltypestring, "zernike")) {
		modeltype = ZERNIKE_LMODEL;
	} else if (!strcmp(modeltypestring, "fourier")) {
		modeltype = FOURIER_LMODEL;
	} else {
		error_exit("readlmodel: unknown model_type\n");
	}

	/* make xitem */
	strptr = (char **) getheaderitemaddress("xname", ipcat);
	if (!strptr) error_exit("readlmodel : missing header item 'xname'\n");
	xname = strptr[0];
	if (strcmp(modeltypestring, "2Dpoly")) {
		dblptr = (double *) getheaderitemaddress("xdim", ipcat);
		if (!dblptr) error_exit("readlmodel : missing header item 'xdim'\n");
		xdim = (int) (dblptr[0]);
	} else {
		xdim = 2;
	}
	xitem = newitem(xname, NUM_TYPE, 1, xdim);

	/* get aitem, aindex, aname */
	if (strcmp(modeltypestring, "2Dpoly")) {
		strptr = (char **) getheaderitemaddress("aname", ipcat);
		if (!strptr) error_exit("readlmodel : missing header item 'aname'\n");
		aname = strptr[0];
		aitem = getobjectitem(aname, ipcat);
		aindex = getobjectitemindex(aname, ipobj);
	} else {
		aindex = 2;
		aitem = (ipcat->itemlist)[aindex];
		aname = aitem->name;
	}
	
	/* and generate the generic model */
	themodel = newlmodel(xitem, aitem);
	themodel->modeltype = modeltype;

	/* copy over the header info */
	copyheaderinfo(themodel->cat, ipcat);

	/* optional xorigin */
	xoriginitem = getheaderitem("xorigin", ipcat);
	if (xoriginitem) {
	        themodel->hasxorigin = 1;
		themodel->xorigin = (double *) xoriginitem->addr;
	}

	/* optional covariance matrix */
	dblptr = (double *) getheaderitemaddress("hascovar", ipcat);
	if (dblptr) {
		if (*dblptr == 1.0) {
			themodel->hascovar = 1.0;
		} else {
			themodel->hascovar = 0.0;
		}
	}


	/* read the cat into a linked list */
	nmodes = 0;
	ipobjbase = NULL;
	while (readobject(ipobj)) {
		ipobj->next = ipobjbase;
		ipobjbase = ipobj;
		ipobj = newobject(ipcat);
		allocobjectcontents(ipobj);
		nmodes++;
	}


	/* set nmodes and allocate lmodel arrays */
	themodel->nmodes = nmodes;
	themodel->a = (void **) calloc(nmodes, sizeof(void *));
	switch (modeltype) {
		case POLYNOMIAL_LMODEL:
			themodel->p = (double **) calloc(nmodes, sizeof(double *));
			themodel->l = (double **) calloc(nmodes, sizeof(double *));
			/* this case is a tad more complicated because we want to be able to read old-style 2Dpoly models */
			if (strcmp(modeltypestring, "2Dpoly")) {
				lindex = getobjectitemindex("l", ipobj);
				pindex = getobjectitemindex("p", ipobj);
			} else {
				lindex = getobjectitemindex("l", ipobj);
				mindex = getobjectitemindex("m", ipobj);
				for (m = 0; m < nmodes; m++) {
					(themodel->p)[m] = (double *) calloc(xdim, sizeof(double));
					(themodel->l)[m] = (double *) calloc(xdim, sizeof(double));
				}
			}
			m = nmodes - 1;
			ipobj = ipobjbase;
			while(ipobj) {
				(themodel->a)[m] = (ipobj->addrlist)[aindex];
				if (strcmp(modeltypestring, "2Dpoly")) {
					(themodel->l)[m] = (ipobj->addrlist)[lindex];
					(themodel->p)[m] = (ipobj->addrlist)[pindex];
				} else {
					ll = ((double **)(ipobj->addrlist))[lindex][0];
					mm = ((double **)(ipobj->addrlist))[mindex][0];
					(themodel->p)[m][0] = ll - mm;
					(themodel->p)[m][1] = mm;
					(themodel->l)[m][0] = ll;
					(themodel->l)[m][1] = mm;
				}
				m--; 
				ipobj = ipobj->next;
			}
			dblptr = (double *) getheaderitemaddress("lmin", ipcat);
			if (dblptr) {
				themodel->lmin = (int) (dblptr[0]);
			}
			dblptr = (double *) getheaderitemaddress("lmax", ipcat);
			if (dblptr) {
				themodel->lmax = (int) (dblptr[0]);
			}
			break;
		case ZERNIKE_LMODEL:
			themodel->n = (double *) calloc(nmodes, sizeof(double));
			themodel->m = (double *) calloc(nmodes, sizeof(double));
			nindex = getobjectitemindex("n", ipobj);
			mindex = getobjectitemindex("m", ipobj);
			m = nmodes - 1;
			ipobj = ipobjbase;
			while(ipobj) {
				(themodel->a)[m] = (ipobj->addrlist)[aindex];
				(themodel->n)[m] = *((double **) (ipobj->addrlist))[nindex];
				(themodel->m)[m] = *((double **) (ipobj->addrlist))[mindex];
				nn = (int) (themodel->n)[m];
				mm = abs((int) (themodel->m)[m]);
				if (mm > nn || 2 * ((nn - mm) / 2) != (nn - mm)) {
					error_exit("readlmodel: illegal m,n indices for zernike polynomial\n");
				}
				m--; 
				ipobj = ipobj->next;
			}
			dblptr = (double *) getheaderitemaddress("nmin", ipcat);
			if (dblptr) {
				themodel->nmin = (int) (dblptr[0]);
			}
			dblptr = (double *) getheaderitemaddress("nmax", ipcat);
			if (dblptr) {
				themodel->nmax = (int) (dblptr[0]);
			}
			makezernikeR(themodel);
			break;
		case FOURIER_LMODEL:
			themodel->k = (double **) calloc(nmodes, sizeof(double *));
			themodel->i = (double *) calloc(nmodes, sizeof(double));
			kindex = getobjectitemindex("k", ipobj);
			iindex = getobjectitemindex("i", ipobj);
			m = nmodes - 1;
			ipobj = ipobjbase;
			while(ipobj) {
				(themodel->a)[m] = (ipobj->addrlist)[aindex];
				(themodel->k)[m] = (ipobj->addrlist)[kindex];
				(themodel->i)[m] = *((double **) (ipobj->addrlist))[iindex];
				m--; 
				ipobj = ipobj->next;
			}
			dblptr = (double *) getheaderitemaddress("kmin", ipcat);
			if (dblptr) {
				themodel->kmin = (int) (dblptr[0]);
			}
			dblptr = (double *) getheaderitemaddress("kmax", ipcat);
			if (dblptr) {
				themodel->kmax = (int) (dblptr[0]);
			}
			dblptr = (double *) getheaderitemaddress("lbox", ipcat);
			if (dblptr) {
				themodel->lbox = dblptr[0];
			} else {
				themodel->lbox = 0.0;
			}
			dblptr = (double *) getheaderitemaddress("boxlimits", ipcat);
			if (dblptr) {
				themodel->boxlimits = (int) (dblptr[0]);
			}
			break;
		default:
			error_exit("readlmodel: unknown model_type\n");
	}
	if (themodel->hascovar) {
		themodel->covar = (double **) calloc(themodel->nmodes, sizeof(double *));
		covarindex = getobjectitemindex("covar", ipobj);
		covaritem = getobjectitem("covar", ipcat);
		if (covaritem->itype != NUM_TYPE || covaritem->ndim != 1 || (covaritem->dim)[0] != themodel->nmodes) {
			error_exit("readlmodel: covar item has wrong type or dimension\n");
		}
		m = nmodes - 1;
		ipobj = ipobjbase;
		while(ipobj) {
			(themodel->covar)[m] = (ipobj->addrlist)[covarindex];
			m--; 
			ipobj = ipobj->next;
		}	
	}
	return(themodel);
}


