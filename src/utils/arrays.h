void	*checkalloc(int nel, int size);
void	allocShortArray(short ***f, int N1, int N2);
void	allocFloatArray(float ***f, int N1, int N2);
void	freeShortArray(short **f, int N1, int N2);
void	freeFloatArray(float **f, int N1, int N2);
void	copyFloatToShort(float **fsrc, short **fdst, int N1, int N2);
void	copyShortToFloat(short **fsrc, float **fdst, int N1, int N2);
float	***alloc3DFloatArray(int N1, int N2, int N3);
void	free3DFloatArray(float ***f, int N1, int N2, int N3);


