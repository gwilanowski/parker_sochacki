function hh_ps=squid_ps
% Function to run the giant squid axon from (Hodgkin and Huxley 1952) implemented
% with the Parker-Sochacki method combined with cubic splines interpolation
% Based on (Stewart and Bair 2009) 
% Written by Grzegorz Wilanowski for Wilanowski, 2012
%clear all; close all; clc; 
warning('off','all'); format long g;
%dt = input('dt = '); %time step size
dt=0.05;  %time step size
%
tol=1e-3; %tolerance
t_end=1000/dt; %end of simulation
%Step current
I=-6.3;   %current is negative
ton=0;
toff=1000;
%Neuron parameters (Hodgkin and Huxley 1952)
cm=1;ena=-115;ek=12;el=-10.613;gna=120;gk=36;gl=0.3;
%Rate constants
v=-150:1:100;
am=0.1*(v+25)./(exp((v+25)./10)-1);
bm=4*exp(v/18);
ah=0.07*exp(v/20);
bh=1./(exp((v+30)/10)+1);
an=0.01*(v+10)./(exp((v+10)./10)-1);
bn=0.125*exp(v/80);
y=[am;bm;ah;bh;an;bn];
%cubic spline interpolation
cs = csapi(v,y);
%Initial values
v=0; %resting potential
am=0.1*(v+25)/(exp((v+25)/10)-1);
bm=4*exp(v/18);
ah=0.07*exp(v/20);
bh=1/(exp((v+30)/10)+1);
an=0.01*(v+10)/(exp((v+10)/10)-1);
bn=0.125*exp(v/80);
%Parameters fo the PS method
fp=zeros(20,1); %floating point parameters
ip=int32(fp); %integer parameters
fp(1)=v;
fp(2)=am/(am+bm);
fp(3)=ah/(ah+bh);
fp(4)=an/(an+bn);
fp(5)=I;
fp(6)=tol;
fp(7)=dt;
fp(8)=ton;
fp(9)=toff;
ip(1)=int32(t_end);
ip(2)=int32(200); %Mclaurin series order
ip(3)=int32(length(cs.breaks)); %spline voltage break points
%Parker-Sochacki algorithm
%cs.coefs', because in MEX columns are exchanged with rows
[v_ps,t_cpu]=ode_ps_hh(fp,ip,cs.coefs',cs.breaks);
t_ps=0:dt:dt*(t_end-1);
%Results
figure;
plot(t_ps,-v_ps);
xlabel('t,s');
ylabel('v, mV');
title('Giant squid axon with Parker Sochacki')
%ISIs
[isi_ps t_ps]=isi(t_ps,v_ps);
ss_ps=mean(isi_ps(0.8*end:end))
disi_ps=ss_ps-isi_ps(1)
mean_ps=mean(isi_ps)
std_ps=std(isi_ps)

fh=figure
plot(t_ps,isi_ps,'*')
xlabel('t,s');
ylabel('isi, s');
title('ISI')
