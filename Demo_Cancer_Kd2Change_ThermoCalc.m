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
Kd2_list=1.1040*10^(-6)*[0, exp(log(10)*linspace(-3.0,3.0,15))] % 20*10^-6 * [0, exp(log(10)*linspace(-1,2,31))]; % weak-binding linker Kd List. Usually recommend to set exponentially-linear spaced values.
N_avogadro=6.02*10^23;

% To change
r_eff = 5.0 % effective radius of weak-binding linker radius, in (nm)

% To change
Cr= 10^2; % Correction ratio, reccommemnd Cr=10^2 to 10^4

V_eff=4/3*pi*(r_eff*10^-9)^3*1000*Cr; % effective volume in (L)
Kd2_eff_list=Kd2_list*N_avogadro*V_eff; % effective Kd2 list

% To change
pA=10^-9; % antibody concentration

% To change
type="randomUniformSphere2D"; % Choose among randomUniformFlat2D randomUniformSphere2D randomQuasiFlat2D randomQuasiSphere2D

% To change 
WperT=2; % number of Weak-binding tether per an antigen

% To change 
isSC=0; % set 1 for considering self-cohesion, 0 for not considering

% To change
L=100; % total area of the surface
density = 1.5; % density of antigen on the surface
Tnum=ceil(L*density)

% So, the total number of antigen (Tnum) becomes L*density

% To change
disp("Parameter setting done")

%% Setting MCMC step Parameters
TestTime=1*2^10;
Project_title = "MMDD_square2D_Tnum98_demo";
IsSave=0; % set 1 to save data, 0 fotherwise
ProbS=zeros(size(Kd2_list,2),TestTime);
 
%% MCMC

parfor i=1:size(Kd2_list,2)
    ProbS(i,:)=par_Metropolis_RS(Project_title,type,L, density,Kd1,Kd2_list(i),Kd2_eff_list(i),pA,TestTime, 10, WperT, isSC)
end
sys_model= Init_AT_System_RS(type,L,density, WperT);

if IsSave
    save("Data\"+Project_title+".mat",'TestTime','Kd2_list','Kd1','V_eff','Kd2_eff_list','pA','sys_model','Tnum','type','ProbS', 'WperT', 'Cr', 'isSC')
end

disp("MCMC done")


%% Drawing figures
cmap=hsv(10);

figure()
dataname=type+" Tnum:"+int2str(Tnum)
semilogx(Kd2_list(2:size(Kd2_list,2)),mean(ProbS(2:size(Kd2_list,2),:),2)/Tnum,'-o','Color',cmap(2,:),'MarkerEdgeColor',cmap(2,:))
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
%% if you want to change log scale to linear scale for x axis, use semilogy() -> it maintain log scale for y-axis only
loglog(Kd2_list(2:size(Kd2_list,2)),pA*(1-mean(ProbS(2:size(Kd2_list,2),:),2)/Tnum)./(mean(ProbS(2:size(Kd2_list,2),:),2)/Tnum),'-o','Color',cmap(2,:),'MarkerEdgeColor',cmap(2,:))
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

%% Drawing with Cohesion
% bounding %
figure();
yyaxis left
data1name="randomUniformSphere2D" % Change dataname
% if you wanna choose different values,
semilogx(Kd2_list(2:size(Kd2_list,2))/k,mean(ProbS(2:size(Kd2_list,2),:),2)/sys_model.Tnum,'-o','Color',cmap(1,:),'MarkerEdgeColor',cmap(1,:))ylabel("log(Kd_eff)");


data2name="SelfCohesion";
data2=[];
for i=2:size(Kd2_list,2)
    data2=[data2 Cohesion(Kd2_list(i)/k,pA,WperT)];
end
yyaxis right
semilogx(Kd2_list(2:size(Kd2_list,2))/k,data2,'-o','Color',cmap(1,:),'MarkerEdgeColor',cmap(2,:))
ylabel("SelfCohesion %")

yyaxis left
data3name="Control"
semilogx([0.5*min(Kd2_list(2:size(Kd2_list,2))),2*max(Kd2_list(2:size(Kd2_list,2)))]/k,[pA/Kd1, pA/Kd1],'Color','blue','LineStyle','--')

title("Kd_2 vs Kd_eff")
xlabel("log(Kd2)")

legend([data1name data3name data2name])

a=[sprintf('%s', datestr(now,'mm-dd-yyyy HH-MM_')) int2str(randi(500))];
if issave
saveas(gcf,"Figure\fig_withCohesion_"+a+".fig")
saveas(gcf,"Figure\fig_withCohesion_"+a+".png")
end



figure();
yyaxis left
data1name="randomUniformSphere2D" % Change dataname
loglog(Kd2_list(2:size(Kd2_list,2)),pA*(1-mean(ProbS(2:size(Kd2_list,2),:),2)/sys_model.Tnum)./(mean(ProbS(2:size(Kd2_list,2),:),2)/sys_model.Tnum),'-o','Color',cmap(1,:),'MarkerEdgeColor',cmap(1,:));hold on;
ylabel("log(Kd_eff)");


data2name="SelfCohesion";
data2=[];
for i=2:size(Kd2_list,2)
    data2=[data2 Cohesion(Kd2_list(i)/k,pA,WperT)];
end
yyaxis right
semilogx(Kd2_list(2:size(Kd2_list,2))/k,data2,'-o','Color',cmap(1,:),'MarkerEdgeColor',cmap(2,:))
ylabel("SelfCohesion %")



yyaxis left
data3name="Control"
semilogx([0.5*min(Kd2_list(2:size(Kd2_list,2))),2*max(Kd2_list(2:size(Kd2_list,2)))],[Kd1, Kd1],'Color','blue','LineStyle','--')

title("Kd_2 vs Kd_eff")
xlabel("log(Kd2)")



legend([data1name data3name data2name])

a=[sprintf('%s', datestr(now,'mm-dd-yyyy HH-MM_')) int2str(randi(500))];

if issave
saveas(gcf,"Figure\fig_withCohesion_"+a+".fig")
saveas(gcf,"Figure\fig_withCohesion_"+a+".png")
end

