/*
 * shmalloc.c
 *
 * routines for allocating and freeing blocks of shared memory
 */


#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <errno.h>
#include "shmalloc.h"

#ifdef DUAL_PROC

int	shmalloc(void **addrptr, int size)
{
	key_t 	key;
	int 	shmflg, shmid;
	void	*addr;

	shmflg = SHM_R | SHM_W;
	key = IPC_PRIVATE;	
	shmid = shmget(key, size, shmflg);
	if (shmid < 0) {
		perror("shmalloc: shmget failed");
		exit(-1);
	}
	addr = shmat(shmid, (char *) 0, shmflg);
	if (addr == (void *) -1) {
		perror("shmalloc: shmat failed");
		exit(-1);
	}
	*addrptr = addr;
	return(shmid);
}

int	shmfree(int shmid)
{
	struct shmid_ds *buf;

	buf = (struct shmid_ds *) calloc(1, sizeof(struct shmid_ds));
	if (shmctl(shmid, IPC_RMID, buf) < 0) {
		perror("shmfree: shmctl failed");
		exit(-1);
	}
	free(buf);
}


int	allocFloatArray_shm(float ***fout, int N1, int N2)
{
	float	*f;
	int	shmid, y;

	shmid = shmalloc((void **) &f, N1 * N2 * sizeof(float));
	*fout = (float **) calloc(N2, sizeof(float *));
	for (y = 0; y < N2; y++) {
		(*fout)[y] = f + N1 * y;
	}
	return (shmid);	
}

#elsif
int	shmalloc(void **addrptr, int size)
{
	error_exit("shmalloc should not be called");
}
int	shmfree(int shmid)
{
	error_exit("shmfree should not be called");
}
int	allocFloatArray_shm(float ***fout, int N1, int N2)
{
	error_exit("shmfree should not be called");
}
#endif
