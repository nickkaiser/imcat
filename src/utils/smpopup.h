/*
 * smpopup.h
 */

void	smpopup(void (*drawfn)(void), void (*cursorfn)(float x, float y, int key), char *geomstring);
void	drawfn(void);
void	cursorfn(float x, float y, int key);

