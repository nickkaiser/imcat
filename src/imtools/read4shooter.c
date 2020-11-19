/*
 * read4shooter.c
 */

#define usage "\n\n\n\
SYNOPSIS\n\
	read4shooter [option...]\n\
\n\
DESCRIPTION\n\
	converts a 4-shooter image to fits format\n\
	Reads and writes from stdin, stdout\n\
\n\
AUTHOR\n\
	Nick Kaiser:  kaiser@cita.utoronto.ca\n\
\n\n\n"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "../imlib/fits.h"
#include "../utils/error.h"
#include "magic.h"
#include "../utils/arrays.h"

#define swap_half(a) ( ((a & 0xff) << 8) | ((unsigned short)(a) >> 8) )
#define N1	816
#define N2	800
#define	NREC	160


short	sheaderrecord();
char	buff[8192];

main(int argc, char *argv[])
{
	int	i, j, rec, comc;
	short	recno, frameno, ccdno, exptime;
	short	uth, utm, uts, rah, ram, ras, rass, decd, decm, decs;
	short	decsign, temperature, filtercode1, filtercode2, base;
	char	objectname[31], date[9], comment[512];
	char	*comv[MAX_COMMENTS];
	short	f[N1];

	if (argc != 1)
		error_exit(usage);
	
	fread(buff, sizeof(char), 8192, stdin);

	recno = sheaderrecord(1);
	frameno = sheaderrecord(3);
	ccdno = sheaderrecord(4);
	exptime = sheaderrecord(5);
	strncpy(objectname, buff + 20, 30);
	objectname[30] = '\0';
	strncpy(date, buff + 70, 8);
	date[8] = '\0';
	uth = sheaderrecord(41);
	utm = sheaderrecord(42);
	uts = sheaderrecord(43);
	rah = sheaderrecord(46);
	ram = sheaderrecord(47);
	ras = sheaderrecord(48);
	rass = sheaderrecord(49);
	decd = sheaderrecord(51);
	decm = sheaderrecord(52);
	decs = sheaderrecord(53);
	decsign = sheaderrecord(55);
	temperature = sheaderrecord(71);
	filtercode1 = sheaderrecord(76);
	filtercode2 = sheaderrecord(77);
	base = sheaderrecord(83);


	comc = 0;
	add_comment(argc, argv, &comc, comv);
	sprintf(comment, "recno = %d", (int) recno);
	add_1comment(comment, &comc, comv);
	sprintf(comment, "frameno = %d", (int) frameno);
	add_1comment(comment, &comc, comv);
	sprintf(comment, "ccdno = %d", (int) ccdno);
	add_1comment(comment, &comc, comv);
	sprintf(comment, "exptime = %d", (int) exptime);
	add_1comment(comment, &comc, comv);
	sprintf(comment, "object: %s", objectname);
	add_1comment(comment, &comc, comv);
	sprintf(comment, "date: %s", date);
	add_1comment(comment, &comc, comv);
	sprintf(comment, "UT = %2d %2d %2d", uth, utm, uts);
	add_1comment(comment, &comc, comv);
	sprintf(comment, "RA = %2d %2d %2d.%d", rah, ram, ras, rass);
	add_1comment(comment, &comc, comv);
	sprintf(comment, "dec = %2d %2d %2d (%1d)", decd, decm, decs, decsign);
	add_1comment(comment, &comc, comv);
	sprintf(comment, "filter = %2d %2d", filtercode1, filtercode2);
	add_1comment(comment, &comc, comv);
	sprintf(comment, "temperature = %d", (int) temperature);
	add_1comment(comment, &comc, comv);
	sprintf(comment, "base = %d", (int) base);
	add_1comment(comment, &comc, comv);

	write_fits_head(N1, N2, comc, comv);
	for (rec = 0; rec < NREC; rec++) {
		if (rec)
			fread(buff, sizeof(char), 8192, stdin);
		else
			for (j = 0; j < N1; j++)
				f[j] = MAGIC;
		for (i = 0; i < 5; i++) {
			for (j = 0; j < N1; j++) {
				f[j] = swap_half(*((short *) (buff + 2 * (N1 * i + j))));
			}
			write_fits_line(f, N1);
		}
	}
	write_fits_tail(N1, N2);
	exit(0);
}


short	sheaderrecord(int i) 
{
	short	result;

	result = *((short *) (buff + 2 * (i - 1)));
	return(swap_half(result));
}
