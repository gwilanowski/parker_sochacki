/*Header file to accompany util.c*/
/*Adapted from Stewart & Bair, 2009*/

#ifndef INC_UTIL_H
#define INC_UTIL_H
typedef unsigned int uint32; /*synonym for unsigned 32 bit integer*/


typedef struct {
	double v;
	double m;
	double h;
	double n;
	double I;
	
	uint32 n_out;
} neuron; 

void first(double **, double **, double *, double **);
void iter(double **, double **, double *, int, double **);
int fabs2(double);

#endif /* INC_UTIL_H */
