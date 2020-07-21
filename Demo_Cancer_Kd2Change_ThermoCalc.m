%%
% Written by Jinyeop Song, 2020/07/20
% This is the demo code for Antibody_ThermoCalc_JY.
% To run the code, follow the description of each sections.

%% Setting System Parameters
clc; clear;

% To change
addpath(genpath('D:\JY_matlab\Antibody_ThermoCalc_JY')) % Add the entire path of Antibody_ThermoCalc_JY

% To change
Kd1=20*10^-9; % antibody - target Kd

% To change
Kd2_list=1.1040*10^(-4)*[0, exp(log(10)*linspace(-2.0,2.0,15))] % 20*10^-6 * [0, exp(log(10)*linspace(-1,2,31))]; % weak-binding linker Kd List. Usually recommend to set exponentially-linear spaced values.
N_avogadro=6.02*10^23;0

% To change
r_eff = 5.0 % effective radius of weak-binding linker radius, in (nm)

V_eff=4/3*pi*(r_eff*10^-9)^3*1000; % effective volume in (L)
Kd2_eff_list=Kd2_list*N_avogadro*V_eff; % effective Kd2 list

% To change
pA=10^-9; % antibody concentration

% To change
type="randomUniformFlat2D"; % Choose among randomUniformFlat2D randomUniformSphere2D randomQuasiFlat2D randomQuasiSphere2D

% To change 
WperT=2;

% To change
L=100; % total area of the surface
density = 1.0; % density of antigen on the surface
% So, the total number of antigen (Tnum) becomes L*density

% To change
disp("Parameter setting done")

%% Setting MCMC step Parameters
TestTime=10*2^10;
Project_title = "MMDD_square2D_Tnum98_demo";
IsSave=0; % set 1 to save data, 0 fotherwise
ProbS=zeros(size(Kd2_list,2),TestTime);
%% to delete
sys=Init_RandomSurface_AT_system(type,L,density, WperT);
sys.Visualize

%%

parfor i=1:size(Kd2_list,2)
    ProbS(i,:)=par_Metropolis_withW(Project_title,type,L, density,Kd1,Kd2_list(i),Kd2_eff_list(i),pA,TestTime, 10, WperT)
end

sys_model=Init_RandomSurface_AT_system(type,Tnum);

if IsSave
    save("Data\"+Project_title+".mat",'TestTime','Kd2_list','Kd1','V_eff','Kd2_eff_list','pA','sys_model','Tnum','type','ProbS')
end

disp("MCMC done")


%% Drawing figures
cmap=hsv(10);

figure()
dataname=type+" Tnum:"+int2str(Tnum)
semilogx(Kd2_list(2:size(Kd2_list,2)),mean(ProbS(2:size(Kd2_list,2),:),2)/sys_model.Tnum,'-o','Color',cmap(2,:),'MarkerEdgeColor',cmap(2,:))

hold on

semilogx([0.5*min(Kd2_list(2:size(Kd2_list,2))),2*max(Kd2_list(2:size(Kd2_list,2)))],[pA/Kd1, pA/Kd1],'Color','blue','LineStyle','--')

title("Kd_2 vs Kd_eff")

xlabel("log(Kd2)")

ylabel("bounding #")



legend([dataname "Control"])

if IsSave
    a=sprintf('%s.', datestr(now,'mm-dd-yyyy HH-MM'));
    saveas(gcf,"Figure\fig_"+a+".fig")
    saveas(gcf,"Figure\fig_"+a+".png")

end



figure()

dataname=type+" Tnum :"+int2str(Tnum)

loglog(Kd2_list(2:size(Kd2_list,2)),pA*(1-mean(ProbS(2:size(Kd2_list,2),:),2)/sys_model.Tnum)./(mean(ProbS(2:size(Kd2_list,2),:),2)/sys_model.Tnum),'-o','Color',cmap(2,:),'MarkerEdgeColor',cmap(2,:))

hold on

semilogx([0.5*min(Kd2_list(2:size(Kd2_list,2))),2*max(Kd2_list(2:size(Kd2_list,2)))],[Kd1, Kd1],'Color','blue','LineStyle','--')

title("Kd_2 vs Kd_eff")

xlabel("log(Kd2)")

ylabel("log(Kd_eff)")



legend([dataname "Control"])

if IsSave

    a=sprintf('%s.', datestr(now,'mm-dd-yyyy HH-MM'));

    saveas(gcf,"Figure\fig_"+a+".fig")

    saveas(gcf,"Figure\fig_"+a+".png")

end



disp("drawing done")



