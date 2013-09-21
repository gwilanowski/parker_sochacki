%Comparsion of ode15s and PS
% clear all
% close all
% clc
load PS
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

