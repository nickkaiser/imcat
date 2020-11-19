/*
 * getop.h - declaration of functions for parsing rpn expressions
 *
 */

#define N_TYPES		5
#define	NUM0_FUNC_TYPE 	0
#define	NUM1_FUNC_TYPE 	1
#define	NUM2_FUNC_TYPE 	2
#define CONSTANT_TYPE	3
#define IM_VALUE_TYPE	4

typedef struct op {
	int	type;
	int	opno;
	double	data;
	void	*addr;
} op;

op	*getop(char **sptr);
char	*getword(char **sptr);
