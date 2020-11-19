/*
 * getop.h - declaration of functions for parsing rpn expressions
 *
 */

item	*getop(char **sptr);
char	*getword(char **sptr);
void	setsourcecathead(cathead *thecathead);
void	setdestcathead(cathead *thecathead);
