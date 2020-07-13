%%
% Written by Jinyeop Song, 2020/07/13 
% This is the demo code for Antibody_ThermoCalc_JY.
% To run the code, follow the description of each sections.

%% Setting System Parameters
clc; clear;

% To change
addpath(genpath('D:\JY_matlab\Antibody_ThermoCalc_JY')) % Add the entire path of Antibody_ThermoCalc_JY

% To change
Kd1=20*10^-9; % antibody - target Kd

% To change
Kd2_list=1.1040*10^(-4)*[0, exp(log(10)*linspace(-3.0,3.0,63)/5)] % 20*10^-6 * [0, exp(log(10)*linspace(-1,2,31))]; % weak-binding linker Kd List. Usually recommend to set exponentially-linear spaced values.

N_avogadro=6.02*10^23;0

% To change
r_eff = 5.0 % effective radius of weak-binding linker radius, in (nm)

V_eff=4/3*pi*(r_eff*10^-9)^3*1000; % effective volume in (L)
Kd2_eff_list=kD2*N_avogadro*V_eff; % effective Kd2 list

% To change
pA=10^-9; % antibody concentration

% To change
type="square2D"; % Choose among "linear1D", "circle1D", "square2D", "traingularSphere2D" ,"cubicSphere2D"

% To change
Tnum=98; % number of antigen number. If Tnum is not allowed by the chosen "type", Tnum is set to nearest allowed value. 

%% Setting MCMC step Parameters
TestTime=50*2^10;
Project_title = "MMDD_square2D_Tnum98_demo";
ProbS=zeros(size(kD2_list,2),TestTime);


%% Run MCMC

if ~isfolder("Data\"+Project_title+"_")
    mkdir("Data\"+Project_title+"_")
end

parfor i=1:size(kD2_list,2)
    ProbS(i,:)=par_Metropolis_withW(Project_title,type,Tnum,kD1,kD2_list(i),kD2_eff_list(i),pA,TestTime, 25)
end

sys_model=Init_AT_System(type,Tnum);
save("Data\"+Project_title+".mat",'TestTime','kD2_list','kD1','V_eff','kD2_eff_list','pA','sys_model','Tnum','type','ProbS')
disp("MCMC done")

%% To be edit, Graph part
Tnum=sys.Tnum;
figure()
semilogx(kD2_list(2:size(kD2_list,2)),mean(ProbS(2:size(kD2_list,2),:),2)/Tnum)
hold on
scatter(kD2_list(2:size(kD2_list,2)),mean(ProbS(2:size(kD2_list,2),:),2)/Tnum,'o','blue')
title(type+" / Kd_2 vs bounding %")
xlabel("log(Kd2)")
ylabel("bounding %")
axis([0.01 1000 0 1])
saveas(gcf,"Figure\"+Project_title+"___kd2vsbounding.png")

figure()
loglog(kD2_list(2:size(kD2_list,2)),(mean(ProbS(2:size(kD2_list,2),:),2)/Tnum)/pA./(1-mean(ProbS(2:size(kD2_list,2),:),2)/Tnum))
hold on
loglog([0.01,1000],[kD1, kD1],'Color','red','LineStyle','--')
hold on
A=max((mean(ProbS(2:size(kD2_list,2),:),2)/Tnum)/pA./(1-mean(ProbS(2:size(kD2_list,2),:),2)/Tnum))
scatter(kD2_list(2:size(kD2_list,2)),(mean(ProbS(2:size(kD2_list,2),:),2)/Tnum)/pA./(1-mean(ProbS(2:size(kD2_list,2),:),2)/Tnum),'o','blue')
hold on
title(type+" / Kd_2 vs Kd_eff")
xlabel("log Kd2")
ylabel("log Kd_eff")
axis([0.01 1000 0 10000+0.1])
saveas(gcf,"Figure\"+Project_title+"___kd2vskdeff.png")

disp("save done")

%% This is function for paraellel computing of MCMC
function ProbS_column = par_Metropolis_withW(Project_title,type,Tnum,kD1,kD2,kD2_eff,pA,TestTime, MCMC_num)

ProbS_column=zeros(1,TestTime);

%disp("start simulation for kD2="+ string(kD2))
sys = Init_AT_System(type,Tnum);
for t=1:TestTime
    if rem(t,2^10)==0
        %disp("Done simulation for kD2="+ string(kD2_list(i))+" & t="+string(t/2^10)+" th 2^10")
    end
    sys = Reinitialize(sys);
    for j=1:MCMC_num
        sys = Metropolis_withW(sys,kD1,kD2,kD2_eff,pA);
    end
    ProbS_column(1,t)=CalculateBindingNum(sys);
end
sys_model=sys;
save("Data\"+Project_title+"_\"+"Kd1_"+string(kD1)+"___Kd2_"+string(kD2)+".mat",'TestTime','kD2','kD1','kD2_eff','pA','sys_model','Tnum','type','ProbS_column')

end