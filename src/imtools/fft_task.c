#include <stdio.h>
#ifdef DUAL_PROC
#include <sys/ipc.h>
#include <sys/shm.h>
#endif
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>

#include "../fftlib/myfft.h"
#include "../utils/error.h"
#include "fft_task.h"

char	*usage = "\n\
NAME\n\
	fft_task --- fft engine for double processor\n\
\n\
SYNOPSYS\n\
	fft_task fr_shmid fk_shmid nx ny mode\n\
\n\
DESCRIPTION\n\
	'fft_task' is invoked to handle 2-D FFTPACK fft's on\n\
	dual processor machines.\n\
	You should not have to invoke it directly.\n\
	It's first two arguments are the id's of shared\n\
	memory segments for the real-space and transformed\n\
	data. nx, ny are dimensions of the fft (must be\n\
	even), and mode can be one of\n\
		0	forward x-transforms\n\
		1	forward y-transforms\n\
		2	inverse x-transforms\n\
		3	inverse y-transforms\n\
\n\
AUTHOR\n\
	Nick Kaiser --- kaiser@ifa.hawaii.edu\n\
\n\n";


main (int argc, char *argv[])
{
	int		fr_shmid, fk_shmid, x, y, nx, ny, mode, shmflg, n;
	int		x1, x2, y1, y2;
	float		*fr;
	complex		*fk;
	float		*r, *work1, *work2;
	pid_t		pid;

	if (argc != 6) {
		fprintf(stderr, usage);
		exit(-1);
	}
	sscanf(argv[1], "%d", &fr_shmid);
	sscanf(argv[2], "%d", &fk_shmid);
	sscanf(argv[3], "%d", &nx);
	sscanf(argv[4], "%d", &ny);
	sscanf(argv[5], "%d", &mode);

	n = ny / 2;

#ifdef DUAL_PROC
	/* attach the data */
	shmflg = SHM_R | SHM_W;
	fr = (float *) shmat(fr_shmid, (char *) 0, shmflg);
	fk = (complex *) shmat(fk_shmid, (char *) 0, shmflg);
	if ((fr == (float *) -1) || (fk == (complex *) -1)) {
		perror("fft_task: shmat failed");
		exit(-1);
	}


	pid = fork();

	if (pid) {
		x1 = 0;
		x2 = nx / 2 - 1;
		y1 = 0;
		y2 = ny / 2 - 1;
	} else {
		x1 = nx / 2;
		x2 = nx - 1;
		y1 = ny / 2;
		y2 = ny - 1;
	}

	switch(mode) {
		case FORWARD_X:
			work2 = (float *) calloc(4 * nx + 15, sizeof(float));
			cffti_(&nx, work2);
			for (y = y1; y <= y2; y++) {
				cfftf_(&nx, fk + nx * y, work2);
			}
			free(work2);
			break;
		case INVERSE_X:
			work2 = (float *) calloc(4 * nx + 15, sizeof(float));
			cffti_(&nx, work2);
			for (y = y1; y <= y2; y++) {
				cfftb_(&nx, fk + nx * y, work2);
			}
			free(work2);
			break;
		case FORWARD_Y:
			r = (float *) calloc(ny, sizeof(float));
			work1 = (float *) calloc(2 * ny + 15, sizeof(float));
			rffti_(&ny, work1);
			for (x = x1; x <= x2; x++) {
				for (y = 0; y < ny; y++) {
					r[y] = fr[y * nx + x];
				}
				rfftf_(&ny, r, work1);
				fk[x].r = r[0];
				fk[x].i = 0.0;
				for (y = 1; y < ny / 2; y++) {
					fk[y * nx + x].r = r[2 * y - 1];
					fk[y * nx + x].i = r[2 * y];
				}
				fk[n * nx + x].r = r[ny - 1];
				fk[n * nx + x].i = 0.0;
			}
			free(work1);
			free(r);
			break;
		case INVERSE_Y:
			r = (float *) calloc(ny, sizeof(float));
			work1 = (float *) calloc(2 * ny + 15, sizeof(float));
			rffti_(&ny, work1);
			for (x = x1; x <= x2; x++) {
				r[0] = fk[x].r;
				for (y = 1; y < n; y++) {
					r[2 * y - 1] = fk[y * nx + x].r;
					r[2 * y] = fk[y * nx + x].i;
				}
				r[ny - 1] = fk[n * nx + x].r;
				rfftb_(&ny, r, work1);
				for (y = 0; y < ny; y++) {
			 		fr[y * nx + x] = r[y] / (nx * ny);
				}
			}
			free(work1);
			free(r);
			break;
		default:
			error_exit("fft_task: bad mode\n");
			break;
	}
	if (pid) {
		wait();
	}
#endif
	exit(0);
}
