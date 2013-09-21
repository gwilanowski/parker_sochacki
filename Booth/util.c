/*Utilities for numerical integration*/
/*Adapted from Stewart & Bair, 2009*/

void first(double **y, double **co, double *fp, double **a){
	double vs,hna,n,mscan,hscan,cas,vd,cad,mcap,mnap,
		mna,mna2,mna3,n2,mscan2,vs2,chis,cas_div,cad_div,
		tauhna,taun,mscan2hscan,vd2,chid;
	double mna3hna,n4,mscan2hscan_vs,chis_vs,chid_vd,
		mcap_vd,iscan,icap,vsp,vdp,I=fp[0];
	double vvs,vsshift=fp[2],vvs2,vvs3,vvd,vdshift=fp[3],vvd2,
		vvd3,hnainf,ninf,mscaninf,hscaninf,mcapinf,mnapinf;
	double co_na,co_kdr,co_scan,co_skca,co_dkca,
		co_cap,co_nap;
	double gc=0.1,ratio=0.1,gna=120,ena=55,gkdr=100,ek=-80,
		gscan=14,eca=80,gskca=3.136,gdkca=.69,gl=0.51,el=-60,                
		f=0.01,alpha1=0.009,alpha2=0.009,kca=2,gcap=0.33,
		gnap=0.2,taumscan=16,tauhscan=160,taumcap=40,
		taumnap=40;   
	vs=y[0][0]; hna=y[1][0]; n=y[2][0]; mscan=y[3][0];
	hscan=y[4][0]; cas=y[5][0]; vd=y[6][0]; cad=y[7][0];
	mcap=y[8][0]; mnap=y[9][0]; mna=y[10][0]; mna2=y[11][0]; 
	mna3=y[12][0]; n2=y[13][0]; mscan2=y[14][0]; vvs2=y[15][0]; 
	chis=y[16][0];cas_div=y[17][0];cad_div=y[18][0];tauhna=y[19][0]; taun=y[20][0];
	mscan2hscan=y[21][0]; vvd2=y[22][0]; chid=y[23][0];
	vvs=vs-vsshift;
	vvs3=vvs2*vvs;
	vvd=vd-vdshift;
	vvd3=vvd2*vvd;
	/*spline interpolation*/
	mna=a[0][0]*vvs3+a[0][1]*vvs2+a[0][2]*vvs+a[0][3];//mna
	tauhna=a[2][0]*vvs3+a[2][1]*vvs2+a[2][2]*vvs+a[2][3];//tauhna
	taun=a[4][0]*vvs3+a[4][1]*vvs2+a[4][2]*vvs+a[4][3];//taun
	mscaninf=a[5][0]*vvs3+a[5][1]*vvs2+a[5][2]*vvs+a[5][3];
	hscaninf=a[6][0]*vvs3+a[6][1]*vvs2+a[6][2]*vvs+a[6][3];
	mcapinf=a[7][0]*vvd3+a[7][1]*vvd2+a[7][2]*vvd+a[7][3];
	mnapinf=a[8][0]*vvd3+a[8][1]*vvd2+a[8][2]*vvd+a[8][3];
	mscan2=mscan*mscan;
    mscan2hscan=mscan2*hscan;		
	mscan2hscan_vs=mscan2hscan*vs;
	iscan=gscan*(mscan2hscan_vs-mscan2hscan*eca);
	y[5][1]=co[5][0]*f*(-alpha1*iscan-kca*y[5][0]);
	cas_div=cas/(cas+0.2);
	n2=n*n;
	n4=n2*n2;
	mna2=mna*mna;
	mna3=mna2*mna;
	mna3hna=mna3*hna;
	co_na=gna*mna3hna;
	co_kdr=gkdr*n4;
	co_scan=gscan*mscan2hscan;
	co_skca=gskca*cas_div;
	chis=-co_na-co_kdr-co_scan-co_skca-gl-gc/ratio;//chis
	chis_vs=chis*vs;
	vsp=co[0][0]*(chis_vs+co_na*ena+co_kdr*ek+
		co_scan*eca+co_skca*ek+gl*el+I+gc*vd/ratio);
	y[0][1]=vsp/1;//vs
	hnainf=a[1][0]*vvs3+a[1][1]*vvs2+a[1][2]*vvs+a[1][3];
	ninf=a[3][0]*vvs3+a[3][1]*vvs2+a[3][2]*vvs+a[3][3];
	y[1][1]=(hnainf-hna)/tauhna;y[1][1]*=co[1][0];
	y[2][1]=(ninf-n)/taun;y[2][1]*=co[2][0];
	y[3][1]=co[3][0]*(mscaninf-mscan)/taumscan;
	y[4][1]=co[4][0]*(hscaninf-hscan)/tauhscan;
	mcap_vd=mcap*vd;
	icap=gcap*(mcap_vd-mcap*eca);
	y[7][1]=co[7][0]*f*(-alpha2*icap-kca*cad);
	cad_div=cad/(cad+0.2);
	co_dkca=gdkca*cad_div;
	co_cap=gcap*mcap;
	co_nap=gnap*mnap;
	//chid=-co_dkca-co_cap-co_nap-gl-gc/(1-ratio);//chid
	chid=-co_dkca-co_cap-co_nap-gl-gc/(1-ratio);//chid
	chid_vd=chid*vd;
	//vdp=co[6][0]*(chid_vd+co_dkca*ek+co_cap*eca+co_nap*ena+
		//gl*el+gc*vs/(1-ratio));
	vdp=co[6][0]*(chid_vd+co_dkca*ek+co_cap*eca+co_nap*ena+
		gl*el+gc*vs/(1-ratio));
	y[6][1]=vdp/1;
	y[8][1]=co[8][0]*(mcapinf-mcap)/taumcap;
	y[9][1]=co[9][0]*(mnapinf-mnap)/taumnap;

}

void iter(double **y, double **co, double *fp, int p, double **a){
	double vs,hna,n,mscan,hscan,cas,vd,cad,mcap,mnap,
		mna,mna2,mna3,n2,mscan2,vs2,chis,cas_div,cad_div,
		tauhna,taun,mscan2hscan,vd2,chid;
	double mna3hna,n4,mscan2hscan_vs,chis_vs,chid_vd,
		mcap_vd,iscan,icap,vsp,vdp,I;
	double vvs,vsshift=fp[2],vvs2,vvs3,vvd,vdshift=fp[3],vvd2,vvd3,
		hnainf,ninf,mscaninf,hscaninf,mcapinf,mnapinf,
		hnainf_next,ninf_next;
	double co_na,co_kdr,co_scan,co_skca,co_dkca,
		co_cap,co_nap;
	double gc=0.1,ratio=0.1,gna=120,ena=55,gkdr=100,ek=-80,
		gscan=14,eca=80,gskca=3.136,gdkca=.69,gl=0.51,el=-60,                
		f=0.01,alpha1=0.009,alpha2=0.009,kca=2,gcap=0.33,
		gnap=0.2,taumscan=16,tauhscan=160,taumcap=40,
		taumnap=40;
	vs=y[0][0]; hna=y[1][0]; n=y[2][0]; mscan=y[3][0];
	hscan=y[4][0]; cas=y[5][0]; vd=y[6][0]; cad=y[7][0];
	mcap=y[8][0]; mnap=y[9][0]; mna=y[10][0]; mna2=y[11][0]; 
	mna3=y[12][0]; n2=y[13][0]; mscan2=y[14][0]; vvs2=y[15][0]; 
	chis=y[16][0];cas_div=y[17][0];cad_div=y[18][0];tauhna=y[19][0]; 
	taun=y[20][0];mscan2hscan=y[21][0]; vvd2=y[22][0]; chid=y[23][0];
	vvs=vs-vsshift;
	vvd=vd-vdshift;
	cauchy_prod(p,y[0],vvs,y[0],vvs,&y[15][p]); //vvs2
	cauchy_prod(p,y[0],vvs,y[15],vvs2,&vvs3); //vvs3
	cauchy_prod(p,y[6],vvd,y[6],vvd,&y[22][p]); //vvd2
	cauchy_prod(p,y[6],vvd,y[22],vvd2,&vvd3); //vvd3
	/*spline interpolation*/
	y[10][p]=a[0][0]*vvs3+a[0][1]*y[15][p]+a[0][2]*y[0][p];//mna
	mscaninf=a[5][0]*vvs3+a[5][1]*y[15][p]+a[5][2]*y[0][p];
	hscaninf=a[6][0]*vvs3+a[6][1]*y[15][p]+a[6][2]*y[0][p];
	mcapinf=a[7][0]*vvd3+a[7][1]*y[22][p]+a[7][2]*y[6][p];
	mnapinf=a[8][0]*vvd3+a[8][1]*y[22][p]+a[8][2]*y[6][p];
	cauchy_prod(p,y[3],mscan,y[3],mscan,&y[14][p]);
	cauchy_prod(p,y[14],mscan2,y[4],hscan,&y[21][p]);//mscan2hscan
	cauchy_prod(p,y[21],mscan2hscan,y[0],vs,&mscan2hscan_vs);
	iscan=gscan*(mscan2hscan_vs-y[21][p]*eca);
	y[5][p+1]=co[5][p]*f*(-alpha1*iscan-kca*y[5][p]);
	series_div(p,y[5][p+1],y[5],cas+0.2,y[17],cas_div);
	cauchy_prod(p,y[2],n,y[2],n,&n2);
	cauchy_prod(p,y[13],n2,y[13],n2,&n4);
	cauchy_prod(p,y[10],mna,y[10],mna,&y[11][p]);//mna2
	cauchy_prod(p,y[10],mna,y[11],mna2,&y[12][p]);//mna3
	cauchy_prod(p,y[12],mna3,y[1],hna,&mna3hna);
	co_na=gna*mna3hna;
	co_kdr=gkdr*n4;
	co_scan=gscan*y[21][p];
	co_skca=gskca*y[17][p];
	y[16][p]=-co_na-co_kdr-co_scan-co_skca;//chis
	cauchy_prod(p,y[16],chis,y[0],vs,&chis_vs);
	vsp=co[0][0]*(chis_vs+co_na*ena+co_kdr*ek+
		co_scan*eca+co_skca*ek+gc*y[6][p]/ratio);
	y[0][p+1]=vsp/(p+1);//vs
	y[19][p]=a[2][0]*vvs3+a[2][1]*y[15][p]+a[2][2]*y[0][p];//tauhna
	y[20][p]=a[4][0]*vvs3+a[4][1]*y[15][p]+a[4][2]*y[0][p];//taun
	cauchy_prod(p+1,y[0],vvs,y[0],vvs,&y[15][p+1]);//vvs2_next
	cauchy_prod(p+1,y[0],vvs,y[15],vvs2,&vvs3);//vvs3_next
	hnainf_next=a[1][0]*vvs3+a[1][1]*y[15][p+1]+a[1][2]*y[0][p+1];
	ninf_next=a[3][0]*vvs3+a[3][1]*y[15][p+1]+a[3][2]*y[0][p+1];
	series_div(p,hnainf_next,y[19],tauhna+1,y[1],hna);y[1][p+1]*=co[1][p];
	series_div(p,ninf_next,y[20],taun+1,y[2],n);y[2][p+1]*=co[2][p];
	y[3][p+1]=co[3][p]*(mscaninf-y[3][p])/taumscan;
	y[4][p+1]=co[4][p]*(hscaninf-y[4][p])/tauhscan;
	cauchy_prod(p,y[8],mcap,y[6],vd,&mcap_vd);
	icap=gcap*(mcap_vd-y[8][p]*eca);
	y[7][p+1]=co[7][p]*f*(-alpha2*icap-kca*y[7][p]);
	series_div(p,y[7][p+1],y[7],cad+0.2,y[18],cad_div);
	co_dkca=gdkca*y[18][p];
	co_cap=gcap*y[8][p];
	co_nap=gnap*y[9][p];
	//y[23][p]=-co_dkca-co_cap-co_nap;//chid
	y[23][p]=-co_dkca-co_cap-co_nap;//chid
	cauchy_prod(p,y[23],chid,y[6],vd,&chid_vd);
	//vdp=co[6][0]*(chid_vd+co_dkca*ek+co_cap*eca+co_nap*ena+
		//gc*y[0][p]/(1-ratio));
	vdp=co[6][0]*(chid_vd+co_dkca*ek+co_cap*eca+co_nap*ena+
		gc*y[0][p]/(1-ratio));
	y[6][p+1]=vdp/(p+1);
	y[8][p+1]=co[8][p]*(mcapinf-y[8][p])/taumcap;
	y[9][p+1]=co[9][p]*(mnapinf-y[9][p])/taumnap;

}

