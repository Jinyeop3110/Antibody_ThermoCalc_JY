%%
% Written by Jinyeop Song, 2020/07/13 
% This is the demo code for Antibody_ThermoCalc_JY.
% To run the code, follow the description of each sections.

%% Setting Parameters
clc; clear;

% To change
addpath(genpath('D:\JY_matlab\Antibody_ThermoCalc_JY')) % Add the entire path of Antibody_ThermoCalc_JY

% To change
kD1=20*10^-9; % antibody - target kD

% To change
kD2=20*10^-5; % weak-binding linker kD

N_avogadro=6.02*10^23;

% To change
r_eff = 5.0 % effective radius of weak-binding linker radius, in (nm)

V_eff=4/3*pi*(r_eff*10^-9)^3*1000; % effective volume in (L)
kD2_eff=kD2*N_avogadro*V_eff; % effective kD2

% To change
pA=10^-9; % antibody concentration

% To change
type="square2D"; % Choose among "linear1D", "circle1D", "square2D", "traingularSphere2D" ,"cubicSphere2D"

% To change
Tnum=98; % number of antigen number. If Tnum is not allowed by the chosen "type", Tnum is set to nearest allowed value. 

%% Run MCMC

sys = Init_AT_System(type,Tnum); % Initiate AT system
Visualize(sys) % Visualization
title("0 th step")

for i=1:5
    sys = Metropolis_withW(sys,kD1,kD2,kD2_eff,pA) % Run the MCMC one step
    sys = Metropolis_withW(sys,kD1,kD2,kD2_eff,pA)
    sys = Metropolis_withW(sys,kD1,kD2,kD2_eff,pA)
    sys = Metropolis_withW(sys,kD1,kD2,kD2_eff,pA)
    sys = Metropolis_withW(sys,kD1,kD2,kD2_eff,pA)
    sys = Metropolis_withW(sys,kD1,kD2,kD2_eff,pA)
    sys = Metropolis_withW(sys,kD1,kD2,kD2_eff,pA)
    sys = Metropolis_withW(sys,kD1,kD2,kD2_eff,pA)
    sys = Metropolis_withW(sys,kD1,kD2,kD2_eff,pA)
    sys = Metropolis_withW(sys,kD1,kD2,kD2_eff,pA)
    sys = Metropolis_withW(sys,kD1,kD2,kD2_eff,pA)

    Visualize(sys)
    title([int2str(i) " * 10 th step"])
end