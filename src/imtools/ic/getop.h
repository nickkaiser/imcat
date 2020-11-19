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

#define N_NUM0_OPS	9
#define	RAND		0
#define X		1
#define Y		2
#define XP		3
#define YP		4
#define IF		5
#define IFQ		6
#define ENTER		7
#define GRAND		8

#define N_NUM1_OPS      22
#define ACOS		0
#define ASIN		1
#define ATAN		2
#define CEIL		3
#define COS		4
#define COSH		5
#define EXP		6
#define FABS		7
#define FLOOR		8
#define LOG		9
#define LOG10		10
#define SIN		11
#define SINH		12
#define SQRT		13
#define TAN		14
#define TANH		15
#define NOT		16
#define BESSJ0		17
#define BESSJ1		18
#define BESSY0		19
#define BESSY1		20
#define SWAP		21



#define N_NUM2_OPS      18
#define	TIMES0		0
#define TIMES1		1
#define PLUS		2
#define DIVIDE		3
#define MINUS		4
#define ATAN2		5
#define POW		6
#define FMOD		7
#define GT		8
#define GE		9
#define LT		10
#define LE		11
#define NE		12		
#define EQ		13		
#define BESSJN		14
#define BESSYN		15
#define MAX		16
#define MIN		17
	

typedef struct op {
	int	type;
	int	opno;
	double	*data;
} op;

void	set_getop_N1(int n1);
op	*getop(char **sptr);
char	*getword(char **sptr);
op	*newop(char *thename);
