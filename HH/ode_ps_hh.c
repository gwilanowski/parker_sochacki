/* =============================================================
 * ode_ps_hh.c Written by Grzegorz Wilanowski, 2012
 * Adapted from (Stewart and Bair 2009)
 *
 * MEX-file for solving the squid giant axon model described 
 * by Hodgkin & Huxley 1952
 * Compile:   mex ode_ps_hh.c util.c integ_util.c
 * Call as:   [v_ps,t_cpu] = ode_ps_hh(fp,ip);
 * ============================================================= */
#include "mex.h"
#include "matrix.h"
#include <math.h> 
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include "util.h"

#define FP     				prhs[0] /*Floating-point parameters*/
#define IP					prhs[1] /*Integer-point parameters*/
#define SP					prhs[2] /*Spline coefs*/
#define BP					prhs[3] /*Spline breaks*/
#define V_PS		      	plhs[0] /*Output voltage*/
#define T_CPU				plhs[1] /*Simulation time*/

#define NV 13 /*Number of variables*/

/*Match voltage v with respective spline breaks bin*/
int sbin(double v,double *bp, int bpsize){
	int k;
	for(k=1;k<bpsize;k++)
		if(v<bp[k]) break;

	return k;
}
 
/*Heaviside function*/
double heav(double x){
    if (x>0) return 1;
	else return 0;
	return -1;
}

int ode_ps_hh(double **yp,double **co,double *yold,double *ynew,neuron *nrn,double *fp,int order_lim,double **a){
	int flag=0,p,nv=4;
	double v,m,h,n,vp2,vp3,n2,n3,n4,m2,m3,m3h,co_k,co_na,gna,gk,gl;
	double psi,ksi,zeta,am,bm,ah,bh,an,bn,chi;
	double eta[NV]; /*Tolerences for v,m,h,n*/
	double tol = fp[1]; /*Tolerance level*/
	double vp; 
	double vshift=fp[2]; /*Beginning of the current splines interval*/

	gna=120;gk=36;gl=0.3;
	
	v = nrn->v; m = nrn->m; h = nrn->h; n = nrn->n; 
	
	/*Get higher power terms*/
	n2=n*n; n4=n2*n2; m2=m*m; m3=m2*m; m3h=m3*h; 
	co_k=gk*n4; co_na=gna*m3h;
	/*Splines interpolation*/
	vp=v-vshift;
	vp2=vp*vp;vp3=vp2*vp;
	am=a[0][0]*vp3+a[0][1]*vp2+a[0][2]*vp+a[0][3];
	bm=a[1][0]*vp3+a[1][1]*vp2+a[1][2]*vp+a[1][3];
	ah=a[2][0]*vp3+a[2][1]*vp2+a[2][2]*vp+a[2][3];
	bh=a[3][0]*vp3+a[3][1]*vp2+a[3][2]*vp+a[3][3];
    an=a[4][0]*vp3+a[4][1]*vp2+a[4][2]*vp+a[4][3];
	bn=a[5][0]*vp3+a[5][1]*vp2+a[5][2]*vp+a[5][3];

	/*Tethered and derived values*/
	psi=-(am+bm);
	ksi=-(ah+bh);
	zeta=-(an+bn);
 	chi = -co_k - co_na - gl;	
 	fp[0]=nrn->I;
 	eta[0]=tol;	eta[1]=tol;	eta[2]=tol;	eta[3]=tol;	
	
	/*Load variables into solver structure*/
 	yp[0][0] = v; yp[1][0] = m; yp[2][0] = h; yp[3][0] = n;
	yp[4][0] = vp2; yp[5][0] = psi;
 	yp[6][0] = ksi; yp[7][0] = zeta; yp[8][0] = m2; yp[9][0] = m3;
 	yp[10][0] = n2; yp[11][0] = chi;

 	/*Run generic PS solver*/
	p = ps_step(yp,co,yold,ynew,fp,eta,first,iter,0,order_lim,nv,a);
 	
	if(p<1){
		flag=1;
	}
	else{		
		nrn->v=ynew[0]; nrn->m=ynew[1]; nrn->h=ynew[2]; nrn->n=ynew[3];
	}
	return flag;
}

/***********************************************************/
void run_sim(double *v_ps,double *t_cpu,double *fp_in,int *ip_in, double *sp, double *bp){
  int t_end=ip_in[0], order_lim=ip_in[1],bpsize=ip_in[2];
  int i,j,k;
  double tol=fp_in[5], dt=fp_in[6];
  double ton=fp_in[7], toff=fp_in[8]; /*Current onset and offset*/ 
  double c0,c1,cp,fp[3];
  double co_v,co_n,co_m,co_h;
  neuron *nrn; 
  int step,p,nv=4,t,flag;  	
  double steps_ps=floor((1.0/dt)+0.5);
  	  
  /*Cell Parameters*/
  double cm=1;
  
  /*Dynamic Data structures for derivs code and generic PS solution*/
  double **a,**co, **yp, *yold, *ynew;
  a = malloc(6*sizeof(double *));
  for(i=0;i<6;i++){
		a[i] = malloc(4*sizeof(double));
  }
  yold = malloc(NV*sizeof(double));  
  ynew = malloc(NV*sizeof(double));
  yp = malloc(NV*sizeof(double *));  
  co = malloc(NV*sizeof(double *));
  for(i=0;i<NV;i++){
	yp[i] = malloc((order_lim+1)*sizeof(double));
	co[i] = malloc((order_lim+1)*sizeof(double));
  }
  nrn = malloc(sizeof(neuron)); 
  	
  fp[1]=tol;
  /*Set variable coefficients*/
  co_v = dt/cm; co_n = dt;	co_m = dt; co_h = dt; /*time rescaling version*/
  co[0][0] = co_v; co[1][0] = co_m; co[2][0] = co_h; co[3][0] = co_n;
	
  /************************************************************/
  /************* Adaptive Parker-Sochacki Method **************/
  /************************************************************/
    
  nrn->v = fp_in[0]; nrn->m = fp_in[1]; nrn->h = fp_in[2];
  nrn->n = fp_in[3];	

  for(p = 1; p < order_lim; p++){ 
	cp = 1.0/(double)(p+1);
	co[0][p] = co[0][0]*cp; co[1][p] = co[1][0]*cp;	co[2][p] = co[2][0]*cp;
	co[3][p] = co[3][0]*cp;	
  }
  c0 = (double)clock();
  /*t is in dt steps, not in seconds!*/
  for(t=0; t<t_end; t++){
  	if ((t*dt>ton)&(t*dt<toff)) 
		/*Step current*/
		nrn->I = fp_in[4];
	else
		nrn->I =0;
	
	v_ps[t] = nrn->v; /*Record voltage*/
	k=sbin(nrn->v,bp,bpsize);
	fp[2]=bp[k-1];
	for(i=0; i<6; i++) /*6 activation/deactivation variables*/
		for(j=0; j<4; j++) /*4 spline interpolation coefs*/
			a[i][j]=sp[4*i+j+24*(k-1)]; /*Spline coefs for bp[k] voltage*/
	
	flag = ode_ps_hh(yp,co,yold,ynew,nrn,fp,order_lim,a); /*PS for t*/
  
  } /*Loop over t*/
  c1 = (double)clock();
  t_cpu[0] = (double)(c1 - c0)/(CLOCKS_PER_SEC);
  mexPrintf("Time = %5.2f. \n",t_cpu[0]); fflush(stdout);
	
	/*Free dynamic arrays*/
  free(nrn); free(yold); free(ynew);
	for(i=0;i<NV;i++){free(yp[i]); free(co[i]);} free(yp); free(co);
    for(i=0;i<6;i++){
		free(a[i]);
	} free(a);
}

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
  double *fp, *v_ps, *t_cpu, dt, *sp,*bp;
  int t_end,*ip;  
  
  /*  Check for proper number of arguments. */
  if (nrhs != 4) mexErrMsgTxt("4 inputs required.");
  if (nlhs != 2) mexErrMsgTxt("2 outputs required.");
	
  /* Get the inputs */ 
  fp = (double *)mxGetData(FP);/*Floating point params*/
  ip = (int *)mxGetData(IP); /*Integer params*/
  t_end=ip[0]; /*Number of time steps, each of size dt*/
  sp = (double *)mxGetData(SP); /*Spline coefs*/
  bp = (double *)mxGetData(BP); /*Spline breaks*/
  
  /* Create output arguments and assign pointers to them */
  V_PS = mxCreateDoubleMatrix(t_end,1,mxREAL); 
  v_ps = (double *)mxGetData(V_PS);
  T_CPU = mxCreateDoubleMatrix(1,1,mxREAL); t_cpu = (double *)mxGetData(T_CPU);

  /* Call the C subroutine. */
  run_sim(v_ps,t_cpu,fp,ip,sp,bp);
}
