/*
 * fitstack.h
 */

int	make_r(void);
int	make_re(void);
int	make_z(void);
int	fitextinctions(void);
int	fittranslations(void);
int	fitdistortions(void);
int	allocatearrays(void);
double	***allocpositionvector(void);
int	printtranslations(char *parfilename);
double	chisquared(void);
int	outputrcat(char *outputcatfilename);
