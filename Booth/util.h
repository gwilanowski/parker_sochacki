/*Header file to accompany tm_util.c*/
/*Adapted from Stewart & Bair, 2009*/

#ifndef INC_UTIL_H
#define INC_UTIL_H
typedef unsigned int uint32; /*synonym for unsigned 32 bit integer*/


typedef struct {
	double vs;
	double hna;
	double n;
	double mscan;
	double hscan;
	double cas;
	double vd;
	double cad;
	double mcap;
	double mnap;
	double I;
	
	uint32 n_out;
} neuron; 

void first(double **, double **, double *, double **);
void iter(double **, double **, double *, int, double **);

#endif /* INC_UTIL_H */
