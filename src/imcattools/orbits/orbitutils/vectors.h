/*
 * vectors.h - declaration of some functions which operate on vectors
 */

void	assign(double *vec, double v0, double v1, double v2);
double	length(double *vec);
void	diff(double *dst, double *src1, double *src2);
void	add(double *dst, double *src1, double *src2);
void	scale(double *vec, double scalefactor);
void	copy(double *dst, double *src);
void	printvec(double *vec, char *name);
void	fprintvec(double *vec, FILE *stream);
double	dot(double *vec1, double *vec2);
void	getperp(double *dst, double *src, double *n);
void	rotx(double *vec, double angle);
void	roty(double *vec, double angle);
void	rotz(double *vec, double angle);
void	cross(double *dst, double *a, double *b);
