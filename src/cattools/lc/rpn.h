/*
 * rpn.h
 */

typedef struct rpnfunc {
	char	*name;
	int	nop;
	item	**oplist;
	item	*result;
}  rpnfunc;

rpnfunc		*newrpnfunction(char *name, char *string);
void		evalrpnfunction(rpnfunc *thefunction);
