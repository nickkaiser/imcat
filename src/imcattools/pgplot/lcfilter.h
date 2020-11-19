/*
 * lcfilter.h
 */

#define NOEDIT_MODE	0
#define REJECT_MODE	1
#define SELECT_MODE	2
#define MASK_MODE	3


void	startfilter(int coordsargcount, char *coordsarg[], int editmode);
int	addfiltercondition(float x1, float x2, float y1, float y2);
void	dofilter(char *tempfilename);
