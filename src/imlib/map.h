/*
 * map.h - declares map(), fastmap(), ultrafastmap(), trianglecount();
 */

#define	ULTRAFAST_MAP_MODE	0
#define FAST_MAP_MODE 		1
#define	TRIANGLE_MAP_MODE	2

 
/*
 * various functions to map a M2 x M1 source image fsource[i][j] onto an N2 x N1 target 
 * image ftarget[i][j] with the mapping 
 *		ftarget(r) = fsource(r + d(r))
 * where the deflection d = (di, dj) is given as a function
 * of target coordinates, is supplied by function deflection()
 * return value of this function is 1/0 for good/MAGIC deflection
 *
 */

/* 
 * map() decomposes target pixel image into triangles
 * slow, but in some sense the best one can do, and can handle
 * inverse mapping - see below.
 */
void	map(float **ftarget, int N1, int N2, float **fsource, int M1, int M2,
			int (*deflection)(float ri, float rj, float *di, float *dj));
long	trianglecount(void);

/*
 * fastmap() uses bi-linear interpolation
 */  
void    fastmap(float **ftarget, int N1, int N2, float **fsource, int M1, int M2,
                        int (*deflection)(float ri, float rj, float *di, float *dj));

/*
 * ultrafastmap() uses nearest pixel
 */  
void    ultrafastmap(float **ftarget, int N1, int N2, float **fsource, int M1, int M2,
                        int (*deflection)(float ri, float rj, float *di, float *dj));


/*
 * 4/19/97: have added these so we can do inverse mapping with map(): The calls
 *	set_triangle_map_mode(INVERSEMAPMODE);
 *	map(fsource, N1, N2, ftarget, M1, M2, def);
 * map the N1 x N2 image fsource onto M1 x M2 image according to
 *	ftarget(x + def(x)) = fsource(x)
 * where the deflection is now given as a function of source coords
 */
#define	FORWARDMAPMODE 1
#define INVERSEMAPMODE 0
void	set_triangle_map_mode(int mapmode);
