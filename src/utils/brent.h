/*
 * brent.h
 *
 * copied from num rec + ansified by Nick Kaiser
 */

float brent(float ax, float bx, float cx, float (*f)(float), 
float tol, float *xmin);
