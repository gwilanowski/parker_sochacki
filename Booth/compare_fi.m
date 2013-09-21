%Comparsion of ode15s and PS
% clear all
% close all
% clc
%
load ODE1
load ODE15s
load PS
load SPLINE
load CNEXP
%
%
scale=0.01;
ton=0;
toff=10000;
tswitch=2500;
%
dt=0.0025; %step size
tbeg=0; %simulation beginning
tfin=10000; %simulation end
t_ps=tbeg:dt:tfin-dt;
%
[f_ode1 t_ode1]=freq(t_ode1,v_ode1);
I_ode1=scale*(t_ode1-ton).*(heav(t_ode1-ton).*heav(toff-t_ode1))+...
    2*scale*(tswitch-t_ode1).*(heav(t_ode1-tswitch).*heav(toff-t_ode1));
[f_ode15s t_ode15s]=freq(t_ode15s,v_ode15s);
I_ode15s=scale*(t_ode15s-ton).*(heav(t_ode15s-ton).*heav(toff-t_ode15s))+...
    2*scale*(tswitch-t_ode15s).*(heav(t_ode15s-tswitch).*heav(toff-t_ode15s));
[f_ps t_ps]=freq(t_ps,v_ps);
I_ps=scale*(t_ps-ton).*(heav(t_ps-ton).*heav(toff-t_ps))+...
    2*scale*(tswitch-t_ps).*(heav(t_ps-tswitch).*heav(toff-t_ps));
[f_spline t_spline]=freq(t_spline,v_spline);
I_spline=scale*(t_spline-ton).*(heav(t_spline-ton).*heav(toff-t_spline))+...
    2*scale*(tswitch-t_spline).*(heav(t_spline-tswitch).*heav(toff-t_spline));
[f_cnexp t_cnexp]=freq(t_cnexp,v_cnexp);
I_cnexp=scale*(t_cnexp-ton).*(heav(t_cnexp-ton).*heav(toff-t_cnexp))+...
    2*scale*(tswitch-t_cnexp).*(heav(t_cnexp-tswitch).*heav(toff-t_cnexp));
%
fh=figure
plot(I_ode1,f_ode1,'k',I_spline,f_spline,'b',I_ps,f_ps,'or',...
I_cnexp,f_cnexp,'m',I_ode15s,f_ode15s,'g')
xlabel('I,\muA');
ylabel('f, Hz');
title('F-I')
fis1=[I_ps' f_ps'];
fis2=[I_cnexp' f_cnexp'];
xlswrite('figures.xls',fis1,'FI1')
xlswrite('figures.xls',fis2,'FI2')
print(fh,'-djpeg',['E:\Grzesiek\Biologia\IBIB\Motoneuron\IVsemestr\sochacki\spline\Article/Figures/'...
    'booth_fi']);
