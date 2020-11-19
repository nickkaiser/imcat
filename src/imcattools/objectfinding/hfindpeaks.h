/*
 * hfindpeaks.h
 */
 
typedef struct peak {
	int		i, j, level;
	float		rf, fs, nu, e1, e2, x[2];
	struct peak	*next;
} peak;

typedef struct line {
	struct peak	*head;
	struct peak	*tail;
	struct line	*next;
	struct line	*prev;
	struct line	*cellmate;
} line;


peak	*getpeaks(int N1, int N2, float **fs, int level, float rf, float sigma);
int	install(peak *thepeak, line *linehead, float dmax);
void	addnewline(line **lineheadptr, peak *thepeak);
int	neighbours(peak *peak1, peak *peak2, float dmax);
void	removedead(line **lineheadptr, int level, float nulimit);
void	disposeof(line *theline, float nulimit);
line	***makecell(line *linehead, float dmax, int ncells);
void	getcoords(peak *thepeak, int *icell, int *jcell, float dmax);
void	freecell(line ***cell, int ncells);
float	thefilterfunction(float ki, float kj);
float	interp(float d, float x1, float x2, float x3);
void	output(peak *bestpk, peak *prevpk, peak *nextpk);