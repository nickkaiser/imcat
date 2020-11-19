/*
 * shmalloc.h
 */

int	shmalloc(void **addrptr, int size);
int	shmfree(int shmid);
int	allocFloatArray_shm(float ***fout, int N1, int N2);
