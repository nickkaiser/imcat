/*
 * lc.h
 */

typedef struct sortobject {
	object	*obj;
	double	sortval;
	struct sortobject *next;
} sortobject;


void	parsename(char *string);
void	parseheaderexpr(char *string);
void	addheadercontents(item *theitem, char *string, cathead *thecathead);
char	*quote(char *string);
void	readtabheader(void);
void	dospecialname(char *string);
void	chopcomments(cathead *thecathead, int mode, int nchop);
void	adddate(void);
int	objcmp(const void *ptr1, const void *ptr2);
