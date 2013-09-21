% Script to run giant squid axon hh model from (Hodgkin and Huxley 1952)
% Written by Grzegorz Wilanowski
clear all; close all; warning('off','all'); format long g;

%Select sims
disp('Starting Hodgkin-Huxley simulation. Select simulation type:')
disp('1 = Giant squid axon');
sim = input('2 = Booth model\n');
if(sim == 1)
    addpath ./HH
    squid_ps
    rmpath('./HH') 
else(sim == 2)
    addpath ./Booth
    booth_ps
    rmpath('./Booth')
end