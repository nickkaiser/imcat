/*
 * map_private.h
 */
 
typedef struct vert {
	float	x, y;
	struct vert	*next, *prev;
} vert;

vert	*makevertex(float x, float y);
void	makering(vert **basevertex, vert **point, int npoints);
void	decompose(vert *basevertex);
void	printpolygon(vert *basevertex);
void	printvertex(vert * thevert);
void	smprintpolygon(vert *basevertex);
void	freepolygon(vert *basevert);
void	dopolygon(vert *basevert);
void	makebreakpoints(vert *v1, vert *v2, vert **V1, vert **V2, int dir, float cut);
void	switchxy(vert *v);
void	link(vert *r1, vert *r2, vert *l1, vert *l2, vert *R1, vert *R2, vert *L1, vert *L2);
void	dotriangle(vert *basevert);
float	trianglearea(vert *basevert);
void	getsourceij(vert *basevert, int *ix, int *iy);
void	addarea(int ix, int iy, float area);
vert	*allocvert(void);
void	freeverts(void);

