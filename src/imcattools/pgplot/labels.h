/*
 * labels.h - stuff for doing labels in contour etc
 */

typedef struct label {
	int	colorindex;
	float	x;
	float	y;
	char	*text;
	struct label *next;
} label;




