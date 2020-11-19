/*
 * linmodel.h
 *
 * functions to fit a stream of data as a linear combination of parameters
 *
 * sample usage:
 *	linmodelinit(np);			allocate space for np-parameter fit
 * 	for (....)
 *		linmodelincrement(f, p);	f = data value: p[0], p[1] ... are parameters
 *	linmodelsolve(F);			fmodel = sum F[i] p[i], free memory
 *
 */
 
 
 void	linmodelinit(int np);
 void	linmodelincrement(double f, double *p);
 void	linmodelsolve(double *F);

