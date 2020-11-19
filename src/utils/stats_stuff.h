typedef struct fstatsrec {
	float 	fmin;
	float 	fmax;
	float 	fmean;
	float 	fmedian;
	float 	fupperquartile;
	float 	flowerquartile;
	float 	fmode;
	float 	sigma;
	long	badpix;
	long	goodpix;
	long	samplesize;
} fstatsrec;

/*
void		do_stats(short **f, int N1, int N2, int margin, statsrec *srec);
*/
void		fdo_stats(float **f, int N1, int N2, int margin, fstatsrec *fsrec);
int		liststats(	float 	*fsample, 
				int	samplesize, 
				float	*median, 
				float	*lquart, 
				float	*uquart, 
				float	*sigma);
int		findmode(float **f, int N1, int N2, float lquart, float uquart, float *mode, float*flquart);



