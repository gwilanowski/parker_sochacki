/* =============================================================
 * ode_ps_booth.c Written by Grzegorz Wilanowski, 2012
 * Adapted from (Stewart and Bair 2009)
 *
 * MEX-file for solving the squid giant axon model described 
 * by Hodgkin & Huxley 1952
 * Compile:   mex ode_ps_booth.c util.c integ_util.c
 * Call as:   [v_ps,t_cpu] = ode_ps_booth(fp,ip);
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

#define NV 24 /*Number of variables*/

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

int ode_ps_booth(double **yp,double **co,double *yold,double *ynew,neuron *nrn,double *fp,int order_lim,double **a){
	int flag=0,p,nv=10;
	int l;
	double vs,hna,n,mscan,hscan,cas,vd,cad,mcap,mnap,
		mna,mna2,mna3,n2,mscan2,vs2,chis,cas_div,cad_div,
		tauhna,taun,mscan2hscan,vd2,chid;
	double mna3hna,n4;
	/*Beginning of the current splines intervals for soma and dendrite*/
	double vsshift=fp[2],vdshift=fp[3];
	double vvs,vvs2,vvs3,vvd,vvd2,vvd3,	
		hnainf,ninf,mscaninf,hscaninf,mcapinf,mnapinf;
	double co_na,co_kdr,co_scan,co_skca,co_dkca,
		co_cap,co_nap;
	double gc=0.1,ratio=0.1,gna=120,ena=55,gkdr=100,ek=-80,
		gscan=14,eca=80,gskca=3.136,gdkca=.69,gl=0.51,el=-60,                
		f=0.01,alpha1=0.009,alpha2=0.009,kca=2,gcap=0.33,
		gnap=0.2;            

	double eta[NV]; /*tolerences for v,m,h,n*/
	double tol = fp[1]; /*tolerance level*/
	
	vs = nrn->vs; hna = nrn->hna; n = nrn->n; mscan = nrn->mscan; 
	hscan = nrn->hscan;cas = nrn->cas;vd = nrn->vd;cad = nrn->cad;
	mcap = nrn->mcap;mnap = nrn->mnap;
	/*Splines interpolation*/
	vvs=vs-vsshift;
	vvs2=vvs*vvs;
	vvs3=vvs2*vvs;
	vvd=vd-vdshift;
	vvd2=vvd*vvd;
	vvd3=vvd2*vvd;
	mna=a[0][0]*vvs3+a[0][1]*vvs2+a[0][2]*vvs+a[0][3];
	hnainf=a[1][0]*vvs3+a[1][1]*vvs2+a[1][2]*vvs+a[1][3];
	tauhna=a[2][0]*vvs3+a[2][1]*vvs2+a[2][2]*vvs+a[2][3];
	ninf=a[3][0]*vvs3+a[3][1]*vvs2+a[3][2]*vvs+a[3][3];
	taun=a[4][0]*vvs3+a[4][1]*vvs2+a[4][2]*vvs+a[4][3];
	mscaninf=a[5][0]*vvs3+a[5][1]*vvs2+a[5][2]*vvs+a[5][3];
	hscaninf=a[6][0]*vvs3+a[6][1]*vvs2+a[6][2]*vvs+a[6][3];
	mcapinf=a[7][0]*vvd3+a[7][1]*vvd2+a[7][2]*vvd+a[7][3];
	mnapinf=a[8][0]*vvd3+a[8][1]*vvd2+a[8][2]*vvd+a[8][3];
	/*Get higher power terms*/
	mna2=mna*mna; mna3=mna2*mna; mna3hna=mna2*hna; 
	n2=n*n; n4=n2*n2; mscan2=mscan*mscan;
	mscan2hscan=mscan2*hscan;
	co_na=gna*mna3hna;
	co_kdr=gkdr*n4;
	co_scan=gscan*mscan2hscan;
	cas_div=cas/(cas+0.2);
	co_skca=gskca*cas_div;
	chis=-co_na-co_kdr-co_scan-co_skca-gl-gc/ratio;
	cad_div=cad/(cad+0.2);
	co_dkca=gdkca*cad_div;
	co_cap=gcap*mcap;
	co_nap=gnap*mnap;
	chid=-co_dkca-co_cap-co_nap-gl-gc/(1-ratio);

 	fp[0]=nrn->I;
	for(l = 0; l < 10; l++){ 
		eta[l] = tol;	
	}
 	/*Load variables into solver structure*/
	yp[0][0]=vs; yp[1][0]=hna; yp[2][0]=n; yp[3][0]=mscan;
	yp[4][0]=hscan; yp[5][0]=cas; yp[6][0]=vd; yp[7][0]=cad;
	yp[8][0]=mcap; yp[9][0]=mnap; yp[10][0]=mna; yp[11][0]=mna2; 
	yp[12][0]=mna3; yp[13][0]=n2; yp[14][0]=mscan2; yp[15][0]=vvs2; 
	yp[16][0]=chis;yp[17][0]=cas_div;yp[18][0]=cad_div;yp[19][0]=tauhna; 
	yp[20][0]=taun;yp[21][0]=mscan2hscan; yp[22][0]=vvd2; yp[23][0]=chid;
	
 	/*Run generic PS solver*/
	p = ps_step(yp,co,yold,ynew,fp,eta,first,iter,0,order_lim,nv,a);
 	
	if(p<1){
		flag=1;
	}
	else{		
		nrn->vs=ynew[0]; nrn->hna=ynew[1]; nrn->n=ynew[2]; 
		nrn->mscan=ynew[3];nrn->hscan=ynew[4]; 
		nrn->cas=ynew[5]; nrn->vd=ynew[6];nrn->cad=ynew[7]; 
		nrn->mcap=ynew[8];nrn->mnap=ynew[9];
	}
	return flag;
}

/***********************************************************/
void run_sim(double *v_ps,double *t_cpu,double *fp_in,int *ip_in, double *sp, double *bp){
  int order_lim=ip_in[1],bpsize=ip_in[2];
  int i,j,k,l,ti;
  double tol=fp_in[11], dt=fp_in[12];
  double ton=fp_in[13], toff=fp_in[14]; /*current onset and offset*/ 
  double c0,c1,cp,fp[4];
  double co_v,co_n,co_m,co_h;
  neuron *nrn; 
  int step,p,flag;  	
  double steps_ps=floor((1.0/dt)+0.5);
  double tswitch=2500, scale=0.01;
  double t,tstart=fp_in[15],tend=fp_in[16];
  int nt=ip_in[0];
    	  
  /*Cell Parameters*/
  double cm=1;
  
  /*Dynamic Data structures for derivs code and generic PS solution*/
  double **a,**co, **yp, *yold, *ynew;
  a = malloc(9*sizeof(double *));
  for(i=0;i<9;i++){
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
  for(l = 0; l < 10; l++){ 
	co[l][0] = dt;	
  }
  co[0][0]=co[0][0]/cm; //co_vs
  co[6][0]=co[6][0]/cm; //co_vd
  
  /************************************************************/
  /************* Adaptive Parker-Sochacki Method **************/
  /************************************************************/
    
  nrn->vs = fp_in[0]; nrn->hna = fp_in[1]; nrn->n = fp_in[2];
  nrn->mscan = fp_in[3]; nrn->hscan = fp_in[4];nrn->cas = fp_in[5];
  nrn->vd = fp_in[6];nrn->cad = fp_in[7];nrn->mcap = fp_in[8];
  nrn->mnap = fp_in[9];

  for(p = 1; p < order_lim; p++){ 
	cp = 1.0/(double)(p+1);
	for(l = 0; l < 10; l++){ 
		co[l][p] = co[l][0]*cp;	
	}
  }

  c0 = (double)clock();
  for(ti=0; ti<nt; ti++){
	  t=tstart+ti*dt;
	  if (fmod((float)ti,1e5)==0){
	  	  mexPrintf("t = %5.2f. \n",t); 
		  mexEvalString("drawnow");
	  }
	  if ((t>=ton)&(t<=toff))
		  nrn->I=scale*(t-ton)*(heav(t-ton)*heav(toff-t))+
		  2*scale*(tswitch-t)*(heav(t-tswitch)*heav(toff-t)); 
		/*Step current*/
		/*nrn->I = fp_in[10];*/
	else
		nrn->I =0;
	
	/*vs spline coefs*/ 
	if(ti>=0){
		v_ps[ti] = nrn->vs;
	}

	k=sbin(nrn->vs,bp,bpsize);
	fp[2]=bp[k-1];
	for(i=0; i<7; i++) /*8 vs activation/deactivation variables*/
		for(j=0; j<4; j++) /*4 spline interpolation coefs*/
			a[i][j]=sp[4*i+j+36*(k-1)]; /*spline coefs for bp[k] voltage*/
	/*vd spline coefs*/ 
	k=sbin(nrn->vd,bp,bpsize);
	fp[3]=bp[k-1];
	for(i=7; i<9; i++) /*2 vd activation/deactivation variables for*/
		for(j=0; j<4; j++) /*4 spline interpolation coefs*/
			a[i][j]=sp[4*i+j+36*(k-1)]; /*spline coefs for bp[k] voltage*/
	

	flag = ode_ps_booth(yp,co,yold,ynew,nrn,fp,order_lim,a); /*PS for t*/
  
  } /*loop over t*/
  c1 = (double)clock();
  t_cpu[0] = (double)(c1 - c0)/(CLOCKS_PER_SEC);
  mexPrintf("Time = %5.2f. \n",t_cpu[0]); fflush(stdout);
	
	/*Free dynamic arrays*/
  free(nrn); free(yold); free(ynew);
	for(i=0;i<NV;i++){free(yp[i]); free(co[i]);} free(yp); free(co);
    for(i=0;i<9;i++){
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
  fp = (double *)mxGetData(FP);/*floating point params*/
  ip = (int *)mxGetData(IP); /*integer params*/
  t_end=ip[0]; /*number of time steps, each of size dt*/
  sp = (double *)mxGetData(SP); /*spline coefs*/
  bp = (double *)mxGetData(BP); /*spline breaks*/
  
  /* Create output arguments and assign pointers to them */
  V_PS = mxCreateDoubleMatrix(t_end,1,mxREAL); 
  v_ps = (double *)mxGetData(V_PS);
  T_CPU = mxCreateDoubleMatrix(1,1,mxREAL); t_cpu = (double *)mxGetData(T_CPU);

  /* Call the C subroutine. */
  run_sim(v_ps,t_cpu,fp,ip,sp,bp);
}
