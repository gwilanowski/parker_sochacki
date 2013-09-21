function hh_ps=booth_ps
% Function to run the Booth model from (Kurian et al. 2010) implemented
% with the Parker-Sochacki method combined with cubic splines interpolation
% Based on (Kurian et al. 2010) and (Stewart and Bair 2009) 
% Written by Grzegorz Wilanowski for Wilanowski, 2012
warning('off','all'); format long g;
%dt = input('dt = '); %time step size
dt=0.005;  %time step size
tol=1e-3; %tolerance
%Soma ramp current
scale=0.01;
ton=0;
toff=10000;
tswitch=2500;
%
tstart=0; %simulation start
tend=10000; %simulation end
%Step current
I=20;   %not used
%Model neuron parameters (Kurian et al. 2010)
gc=0.1;
p=0.1;
c=1;
gna=120;
ena=55;
thetamna=-35;           
kmna=-7.8;
thetahna=-55;
khna=7;
gkdr=100;
ek=-80;
thetan=-28;
kn=-15;
gscan=14;
eca=80;
thetamsca=-30;
kmscan=-5;
taumscan=16;            
thetahsca=-45;
khscan=5;
tauhscan=160;           
gskca=3.136;
gdkca=.69;
gl=0.51;
el=-60;                
f=0.01;
alpha1=0.009;
alpha2=0.009;
kca=2;
gcap=0.33;             
thetamcap=-40;
kmcap=-7;
taumcap=40;
gnap=0.2;               
thetamnap=-25;
kmnap=-4;
taumnap=40;
%Set up vector of variables
soma_var=6;             % First 6 entries are for soma variables
dend_var=4;             % # of variables in dendrite compartment
sy=soma_var + dend_var; % # of equations
%Variables to be interpolated with spline
v=-150:1:100;
mna=1./(1+exp((v-thetamna)/kmna));
hnainf=1./(1+exp((v-thetahna)/khna));
tauhna=120./(exp((v+50)/15)+exp(-(v+50)/16));
ninf=1./(1+exp((v-thetan)/kn));
taun=28./(exp((v+40)/40)+exp(-(v+40)/50));
mscaninf=1./(1+exp((v-thetamsca)/kmscan));
hscaninf=1./(1+exp((v-thetahsca)/khscan));
mcapinf=1./(1+exp((v-thetamcap)/kmcap));
mnapinf=1./(1+exp((v-thetamnap)/kmnap));
%Cubic spline interpolation
y=[mna;hnainf;tauhna;ninf;taun;mscaninf;hscaninf;mcapinf;mnapinf];
cs = csapi(v,y);
%Parameters fo the PS method
fp=zeros(20,1); %floating point parameters
ip=int32(fp); %integer parameters
fp(11)=I;
fp(12)=tol;
fp(13)=dt;
fp(14)=ton;
fp(15)=toff;
ip(2)=int32(200); %Mclaurin series order
ip(3)=int32(length(cs.breaks)); %spline voltage break points
%
FP=load('FP.mat');
fp(1:10)=FP.fp;
fp(16)=tstart;
fp(17)=tend;
nt=floor((tend-tstart)/dt+0.5); %no of time steps
ip(1)=int32(nt); 
%cs.coefs', because in MEX columns are exchanged with rows
[v_ps,t_cpu]=ode_ps_booth(fp,ip,cs.coefs',cs.breaks);
%Plot time response
figure;
t_ps=tstart:dt:tend-dt;
plot(t_ps,v_ps);
xlabel('t,s');
ylabel('v, mV');
title('Booth model with Parker Sochacki')
%Plot f-I curve
figure
[f_ps t_ps]=freq(t_ps,v_ps);
I_ps=scale*(t_ps-ton).*(heav(t_ps-ton).*heav(toff-t_ps))+...
    2*scale*(tswitch-t_ps).*(heav(t_ps-tswitch).*heav(toff-t_ps));
plot(I_ps,f_ps,'*')
xlabel('I,\muA');
ylabel('f, Hz');
title('f-I')
%
load gong;
player = audioplayer(y, Fs);
play(player);