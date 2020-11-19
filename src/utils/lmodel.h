/* 
 * lmodel.h -- declare lmodel structure and related functions
 */

typedef struct lmodel {
	/* generic lmodel elements */
	int	modeltype;
	cathead	*cat;
	item	*xitem, *aitem;
	int 	xdim;
	char	*xname, *aname;
	int 	nmodes, asize, j;
	int	hasxorigin;
	double	*xorigin;
	void	**a;
	double	**aflat;
	double	hascovar;
	double	**covar;
	/* polynomial stuff */
	int	lmin, lmax;
	int	*pp;
	double	**p, **l;
	/* zernike stuff */
	int	nmin, nmax;
	double	*n, *m, **R;
	/* fourier stuff */
	int	kmin, kmax, boxlimits;
	int	*kk;
	double	**k, *i, lbox;
} lmodel;

#define POLYNOMIAL_LMODEL 	0
#define ZERNIKE_LMODEL 		1
#define	FOURIER_LMODEL		2



/* lmodel.c */

/* newlmodel only creates modeltype independent parts of lmodel structure */
lmodel *newlmodel(item *xitem, item *aitem);



/* create a full polymodel l-model structure */
lmodel *newpolylmodel(item *xitem, item *aitem, int lmin, int lmax);

/* recursive function to loop over all polynomial modes of order l */
int	polyloop(lmodel *themodel, int l, int i);



/* create a full zernike l-model structure */
lmodel *newzernikelmodel(item *xitem, item *aitem, int nmin, int nmax);

/* recursive function to loop over all zernike modes of order n */
int	zernikeloop(lmodel *themodel, int n, int i);

/* generate the R array */
int	makezernikeR(lmodel *themodel);



/* create a full fourier l-model structure */
lmodel *newfourierlmodel(item *xitem, item *aitem, int kmin, int kmax, int boxlimits, double lbox);

/* recursive function to loop over all fourier modes */
int	fourierloop(lmodel *themodel, int i);



/* the mode functions */
double	lmodelfunc(lmodel *themodel, int m, double *x);

/* functions to convert from a[m][][]... to aflat[j][m] and vice versa */
int	flatten_a(lmodel *themodel);
int	rflatten_a(lmodel *themodel, void *a, int m, int level);
int	unflatten_a(lmodel *themodel);
int	runflatten_a(lmodel *themodel, void *a, int m, int level);

/* recursive function to add fac times asrc to adst */
int	addtomatrix(void *adst, void *asrc, double fac, int ndim, int *dim, int level);

/* recursive function to zero a matrix */
int	zeromatrix(void *a, int ndim, int *dim, int level);

/* factorial */
int	bang(int n);



/* lmodelio.c */

/* write out a lmodel structure as lc-format cat */
int	writelmodel(lmodel *themodel, FILE *opstream);

/* read a lmodel structure */
lmodel	*readlmodel(FILE *ipstream);




/* lmodelcalculus.c */

/* create derivative of an lmodel */
lmodel	*difflmodel(lmodel *iplmodel, int integratemode);
