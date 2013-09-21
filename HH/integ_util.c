/*Integration routine for the Parker-Sochacki method*/

/*Parker-Sochacki stepper - Solve and update in one function*/
int ps_step(double **y,double **co,double *y1,double *ynew,double *fp,
  double eta[],void (*first)(double **,double **,double *), 
	void (*iter)(double **,double **,double *,int),int stop,int ps_limit, int nv, double **a){
  int i,p; 
  first(y,co,fp,a); /*Calculate first order terms*/
  for(i=0;i<nv;i++)y1[i]=y[i][0]+y[i][1];
  for(p=1;p<(ps_limit-1);p++){/*Iterations*/
    iter(y,co,fp,p,a);
	for(i=0;i<nv;i++)ynew[i]=y1[i]+y[i][p+1]; /*Update*/
    if((y[0][p+1]>10.0)||(y[0][p+1]<-10.0)){p=-1;break;} /*Check for divergence*/
	for(i=0;i<nv;i++){ /*Check for tolerance*/
		if(((ynew[i]-y1[i])<eta[i])||((ynew[i]-y1[i])>-eta[i]))
			break; /*if tolerance is satisfied*/
	}
	for(i=0;i<nv;i++)y1[i]=ynew[i];
  } p++;
  return p;
}

/*Cauchy product calculation from Stewart & Bair, 2009*/
void cauchy_prod(int p,double *a,double a0,double *b,double b0,double *c){
	/*c is the pth term of the Cauchy product of a and b with zeroth order 
	terms a0, b0 allowing shifted products (e.g (a-1)*(b+1))*/
 	int i;
 	*c = a0*b[p] + b0*a[p];
	for(i = 1; i < p; i++){*c += a[i]*b[p-i];}
}