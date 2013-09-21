/*Utilities for numerical integration*/
/*Adapted from Stewart & Bair, 2009*/

void first(double **y, double **co, double *fp, double **a){
	double I,v,m,h,n,vp2,psi,ksi,zeta,m2,m3,n2,chi,vp3,m3h,n4;
	double am,bm,ah,bh,an,bn,co_na,co_k,gna,gk,gl,ena,ek,el;
	double v1;double vp;double vshift=fp[2];
	ena=-115;ek=12;el=-10.613;gna=120;gk=36;gl=0.3;
	I=fp[0];
	v=y[0][0]; m=y[1][0]; h=y[2][0]; n=y[3][0];
	vp2=y[4][0]; psi=y[5][0]; ksi=y[6][0]; zeta=y[7][0];
	m2=y[8][0]; m3=y[9][0]; n2=y[10][0]; chi=y[11][0];
	vp=v-vshift;
	vp3=vp2*vp;
	/*spline interpolation*/
	am=a[0][0]*vp3+a[0][1]*vp2+a[0][2]*vp+a[0][3];
	bm=a[1][0]*vp3+a[1][1]*vp2+a[1][2]*vp+a[1][3];
	ah=a[2][0]*vp3+a[2][1]*vp2+a[2][2]*vp+a[2][3];
	bh=a[3][0]*vp3+a[3][1]*vp2+a[3][2]*vp+a[3][3];
    an=a[4][0]*vp3+a[4][1]*vp2+a[4][2]*vp+a[4][3];
	bn=a[5][0]*vp3+a[5][1]*vp2+a[5][2]*vp+a[5][3];

	m3h=m3*h;
    n4=n2*n2;
	co_na=gna*m3h;
	co_k=gk*n4;
	v1=co[0][0]*(v*chi+co_na*ena+co_k*ek+gl*el+I); //v'
	y[0][1]=v1;
	y[1][1]=co[1][0]*(m*psi+am);
	y[2][1]=co[2][0]*(h*ksi+ah);
	y[3][1]=co[3][0]*(n*zeta+an);

}

void iter(double **y, double **co, double *fp, int p, double **a){
	double v,m,h,n,vp2,psi,ksi,zeta,m2,m3,n2,chi,vp3,m3h,n4;
	double am,bm,ah,bh,an,bn,co_na,co_k,gna,gk,gl,ena,ek;
	double m_psi,h_ksi,n_zeta,v_chi;
	double vp;double vshift=fp[2];
	ena=-115;ek=12;gna=120;gk=36;gl=0.3;
		
	v=y[0][0]; m=y[1][0]; h=y[2][0]; n=y[3][0];
	vp2=y[4][0]; psi=y[5][0]; ksi=y[6][0]; zeta=y[7][0];
	m2=y[8][0]; m3=y[9][0]; n2=y[10][0]; chi=y[11][0];
	vp=v-vshift;
	cauchy_prod(p,y[0],vp,y[0],vp,&y[4][p]); //vp2
	cauchy_prod(p,y[0],vp,y[4],vp2,&vp3); //vp3
	
	/*spline interpolation*/
	am=a[0][0]*vp3+a[0][1]*y[4][p]+a[0][2]*y[0][p];
	bm=a[1][0]*vp3+a[1][1]*y[4][p]+a[1][2]*y[0][p];
	ah=a[2][0]*vp3+a[2][1]*y[4][p]+a[2][2]*y[0][p];
	bh=a[3][0]*vp3+a[3][1]*y[4][p]+a[3][2]*y[0][p];
    an=a[4][0]*vp3+a[4][1]*y[4][p]+a[4][2]*y[0][p];
	bn=a[5][0]*vp3+a[5][1]*y[4][p]+a[5][2]*y[0][p];
	
	
	y[5][p]=-(am+bm); //psi
	y[6][p]=-(ah+bh); //ksi
	y[7][p]=-(an+bn); //zeta
	cauchy_prod(p,y[5],psi,y[1],m,&m_psi); //m*psi
	cauchy_prod(p,y[6],ksi,y[2],h,&h_ksi); //h*ksi
	cauchy_prod(p,y[7],zeta,y[3],n,&n_zeta); //n*zeta
	cauchy_prod(p,y[1],m,y[1],m,&y[8][p]); //m2
	cauchy_prod(p,y[1],m,y[8],m*m,&y[9][p]); //m3
	cauchy_prod(p,y[2],h,y[9],m*m*m,&m3h); //m3h
	cauchy_prod(p,y[3],n,y[3],n,&y[10][p]); //n2
	cauchy_prod(p,y[10],n*n,y[10],n*n,&n4); //n4
	co_na=gna*m3h;
	co_k=gk*n4;
	y[11][p]=-co_na-co_k; //chi
	cauchy_prod(p,y[11],chi,y[0],v,&v_chi); //v*chi
	y[12][p]=co[0][0]*(v_chi+co_na*ena+co_k*ek); //v'

	y[0][p+1]=y[12][p]/(p+1);
	y[1][p+1]=co[1][p]*(m_psi+am);
	y[2][p+1]=co[2][p]*(h_ksi+ah);
	y[3][p+1]=co[3][p]*(n_zeta+an);

}

