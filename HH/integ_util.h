#ifndef INC_INTEG_UTIL_H
#define INC_INTEG_UTIL_H

extern int ps_step(double **,double **,double *,double *,double *,
  double *,void (*)(double **,double **,double *), 
	void (*)(double **,double **,double *,int), int, int, int);

extern void cauchy_prod(int, double *, double, double *, double, double *);


#endif /* INC_INTEG_UTIL_H */
