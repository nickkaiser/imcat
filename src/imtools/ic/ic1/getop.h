/*
 * getop.h - declaration of functions for parsing rpn expressions
 *
 */

#define	NUM0_FUNC_TYPE 	4
#define	NUM1_FUNC_TYPE 	0
#define	NUM2_FUNC_TYPE 	1
#define CONSTANT_TYPE	2
#define IM_VALUE_TYPE	3

typedef struct op {
	int	type;
	int	opno;
	double	data;
} op;

op	*getop(char **sptr);
char	*getword(char **sptr);
