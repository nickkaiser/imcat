/*
 * fft_sgi.h
 */
typedef struct transform {
	float	***data;
	float	**speq;
} transform;
	
typedef	transform *fft_type;
