/*
 * getpsf.c
 */


#define	usage "\n\n\n\
SYNOPSIS\n\
		smselect xname yname [options...] [cat]\n\
			-X x1 x2	# specify range for X-coord\n\
			-Y y1 y2	# specify range for Y-coord\n\
			-v		# select objects outside the box\n\
			-d		# just display the plot\n\
			-m		# return average x and y\n\
			-V		# verbose mode\n\
\n\
DESCRIPTION\n\
		\"smselect\" interactive 2-d catalogue editor.\n\
		Invokes sm to popup a window with the objects plotted.\n\
		User then generates a box using cursor and key-'b'\n\
		and when this is OK selects objects lying within the\n\
		box by hitting key-'x'.\n\
		Use -v to select objects which don't lie in the box.\n\
		Reads an lc-format catalogue 'cat' and writes to stdout.\n\
		You need to quote fancy names (e.g. 'x[0]').\n\
\n\n\n"		


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../imlib/fits.h"
#include "../utils/error.h"
#include "../utils/smpopup.h"

#define       MAX(x,y) (((x) > (y)) ? (x) : (y))
#define       MIN(x,y) (((x) < (y)) ? (x) : (y))

float		*smx, *smy, smx1, smy1, smx2, smy2;
float		smxlim1, smylim1, smxlim2, smylim2;
char		*smxname, *smyname;
int		smN;
void		drawselectionbox();
void		getlimits(int N, float *x, float *xmin, float *xmax);

main(int argc, char *argv[])	
{
	char		line[64], systemstring[1024], *catfilename, tempfile1[128], tempfile2[128];
	int		arg, exclude, needxlimits, needylimits, displayonly, copycat, doaverage, verbose;
	float		xmin, xmax, ymin, ymax, xbar, ybar;
	int		i, nobjects, pid;
	FILE		*tempf, *pipe;

	/* defaults */
	exclude = 0;
	needxlimits = 1;
	needylimits = 1;
	displayonly = 0;
	doaverage = 0;
	verbose = 0;
	
	
	/* parse args */
	if (argc < 3)
		error_exit(usage);
	smxname = argv[1];
	smyname = argv[2];
	arg = 3;	
	while (arg < argc) {
		if (argv[arg][0] != '-') {
			break;
		}
		switch (argv[arg++][1]) {
			case 'X':
				sscanf(argv[arg++], "%f", &smxlim1);
				sscanf(argv[arg++], "%f", &smxlim2);
				needxlimits = 0;
				break;
			case 'Y':
				sscanf(argv[arg++], "%f", &smylim1);
				sscanf(argv[arg++], "%f", &smylim2);
				needylimits = 0;
				break;
			case 'v':
				exclude = 1;
				break;
			case 'd':
				displayonly = 1;
				verbose = 0;
				break;
			case 'V':
				verbose = 1;
				break;
			case 'm':
				doaverage = 1;
				verbose = 0;
				break;
			default:
				error_exit(usage);
				break;
		}
	}

	/* make temporary file names */
	pid = getpid();
	sprintf(tempfile1, "%dA.tmp", pid);
	sprintf(tempfile2, "%dB.tmp", pid);


	if (arg == argc) {
		sprintf(systemstring, "cat > %s", tempfile1);
		pipe = popen(systemstring, "r");
		pclose(pipe);
		catfilename = tempfile1;
	} else {
		catfilename = argv[arg];
	}

	sprintf(systemstring, "lc -o 'x = %%%s' 'y = %%%s' < %s > %s", smxname, smyname, catfilename, tempfile2);
	if (system(systemstring))
		error_exit("smselect: lc command failed\n");
	
	/* count the objects */
	tempf = fopen(tempfile2, "r");
	nobjects = 0;
	while (fgets(line, 63, tempf))
		nobjects++;
	fclose(tempf);
	
	/* allocate and fill the x, y arrays */
	smx = (float *) calloc(nobjects, sizeof(float));
	smy = (float *) calloc(nobjects, sizeof(float));
	smN = nobjects;

	/* read the objects */
	tempf = fopen(tempfile2, "r");
	nobjects = 0;
	while (fgets(line, 63, tempf))
		sscanf(line, "%f %f", smx + nobjects, smy + nobjects++);
	fclose(tempf);

	if (needxlimits)
		getlimits(smN, smx, &smxlim1, &smxlim2);
	if (needylimits)
		getlimits(smN, smy, &smylim1, &smylim2);


	/* make x y scatter plot & get user to select stars */
	if (verbose) {
		fprintf(stderr, "#\n# use cursor + key 'b' until you get a nice box around stars\n");
		fprintf(stderr, "# then hit 'x' to exit from cursor input mode\n#\n");
	}
	smpopup(drawfn, cursorfn, "-g 400x400+100+100");
	xmin = MIN(smx1, smx2);
	xmax = MAX(smx1, smx2);
	ymin = MIN(smy1, smy2);
	ymax = MAX(smy1, smy2);
	if (verbose) {
	fprintf(stderr, "# selected range\n#\tx = %10.3e to %10.3e\n#\ty = %10.3e to %10.3e\n",
		xmin, xmax, ymin, ymax);
	}

	if (doaverage) {
		xbar = ybar = 0.0;
		nobjects = 0;
		for (i = 0; i < smN; i++) {
			if (smx[i] >= xmin && smx[i] <= xmax && smy[i] >= ymin && smy[i] <= ymax) {
				if (!exclude) {
					xbar += smx[i];
					ybar += smy[i];
					nobjects ++;
				}
			} else {
				if (exclude) {
					xbar += smx[i];
					ybar += smy[i];
					nobjects ++;
				}
			}
		}
		if (nobjects) {
			xbar /= nobjects;
			ybar /= nobjects;
		}
		sprintf(systemstring, "rm -f %s %s", tempfile1, tempfile2);
		system(systemstring);
		fprintf(stdout, "%f %f\n", xbar, ybar);
		exit(0);
	}

	if (!displayonly) {
		if (exclude) {
			sprintf(systemstring, "lc -b -i '%%%s %f > %%%s %f < and %%%s %f > and %%%s %f < and !' < %s ", 
				smxname, xmin, smxname, xmax, smyname, ymin, smyname, ymax, catfilename);
		} else {
			sprintf(systemstring, "lc -b -i '%%%s %f > %%%s %f < and %%%s %f > and %%%s %f < and' < %s ", 
				smxname, xmin, smxname, xmax, smyname, ymin, smyname, ymax, catfilename);
		}
		if (system(systemstring))
			error_exit("smselect: lc command failed\n");
	}
	sprintf(systemstring, "rm -f %s %s", tempfile1, tempfile2);
	system(systemstring);
	exit(0);
}


void	drawfn(void)
{
	float	dottype = 11.0, crosstype = 41.0, startype = 62.0;

	sm_limits(smxlim1, smxlim2, smylim1, smylim2);
	sm_box(1, 2, 0, 0);
	sm_ptype(&dottype, 1);
	sm_points(smx, smy, smN);
	sm_xlabel(smxname);
	sm_ylabel(smyname);
}




void	cursorfn(float x, float y, int key)
{
	static	int first = 1;

	if (key == 'b') {
		if (first) {
			smx1 = x;
			smy1 = y;
			first = 0;
		} else {
			smx2 = x;
			smy2 = y;
			sm_erase();
			drawfn();
			drawselectionbox();
			first = 1;
		}
	}
}





void	drawselectionbox(void)
{
	sm_relocate((float) smx1, (float) smy1);
	sm_draw((float) smx1, (float) smy2);
	sm_draw((float) smx2, (float) smy2);
	sm_draw((float) smx2, (float) smy1);
	sm_draw((float) smx1, (float) smy1);
	sm_gflush();
}



#define MARGIN 0.1

void		getlimits(int N, float *x, float *xmin, float *xmax)
{
	int	i;
	float	d;

	*xmin = *xmax = x[0];
	for (i = 1; i < N; i++) {
		*xmin = (x[i] < *xmin ? x[i] : *xmin);
		*xmax = (x[i] > *xmax ? x[i] : *xmax);
	}
	d = *xmax - *xmin;
	*xmax += MARGIN * d;
	*xmin -= MARGIN * d;	
}
#undef MARGIN