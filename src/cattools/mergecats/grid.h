/*
 * grid.h
 */

void	setgridsize(double *xmin, double *xmax, double  d);
void	allocgrid(object ****gridptr);
int	getgridcoords(double *x, int *ix, int *iy);
void	getneighbours(double *x, object ***grid, int *nneighbours, object **neighbourobj);
void	getobjects(object ***grid, int *nobjs, object **theobjlist);


#define MAX_NEIGHBOURS 100000
