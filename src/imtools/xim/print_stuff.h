/* public use: */
void	print_im(int N1, int N2, int Ncolors, double (*gray)(), char *caption);	
void	setpagesize(int pagewidth, int pageheight, int pagemargin, int dosetpagedevice);
	


/* private use: */
void	print_caption(char *caption);
void	psDrawChar(char aChar);
void	ps(char *string);
void	set_print_opf(FILE *thefile);







