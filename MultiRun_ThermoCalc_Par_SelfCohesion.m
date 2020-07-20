%%
clc; clear;
%% Setting System Param
Kd1=20*10^-9; % antibody - target
Kd2_list=2.0*10^(-6)*[0, exp(log(10)*linspace(-4.0,4.0,127))]; % 20*10^-6 * [0, exp(log(10)*linspace(-1,2,31))]; % weak-binding
N_avogadro=6.02*10^23;
k=10^2;
V_eff=4/3*pi*(50*10^-10)^3*1000*k;
Kd2_eff_list=Kd2_list*N_avogadro*V_eff;
pA=10^-9;
type="cubicSphere2D"; % "linear1D", "circle1D", "square2D", "traingularSphere2D" ,"cubicSphere2D"
Tnum=98;
%% Setting Exp param
TestTime=20*2^10;
Project_title = "0720_squareSphere2D_Tnum98_AfterCorrection";
ProbS=zeros(size(Kd2_list,2),TestTime);

%% Initiate AT system

if ~isfolder("Data\"+Project_title+"_")
    mkdir("Data\"+Project_title+"_")
end

parfor i=1:size(Kd2_list,2)
    ProbS(i,:)=par_Metropolis_withW_temp(Project_title,type,Tnum,Kd1,Kd2_list(i),Kd2_eff_list(i),pA,TestTime, 20)
end

sys_model=Init_AT_System(type,Tnum);
save("Data\"+Project_title+".mat",'TestTime','Kd2_list','Kd1','V_eff','Kd2_eff_list','pA','sys_model','Tnum','type','ProbS')
disp("MCMC done")
%% 
Tnum=sys.Tnum;
figure()
semilogx(Kd2_list(2:size(Kd2_list,2)),mean(ProbS(2:size(Kd2_list,2),:),2)/Tnum)
hold on
scatter(Kd2_list(2:size(Kd2_list,2)),mean(ProbS(2:size(Kd2_list,2),:),2)/Tnum,'o','blue')
title(type+" / Kd_2 vs bounding %")
xlabel("log(Kd2)")
ylabel("bounding %")
axis([0.01 1000 0 1])
saveas(gcf,"Figure\"+Project_title+"___kd2vsbounding.png")

figure()
loglog(Kd2_list(2:size(Kd2_list,2)),(mean(ProbS(2:size(Kd2_list,2),:),2)/Tnum)/pA./(1-mean(ProbS(2:size(Kd2_list,2),:),2)/Tnum))
hold on
loglog([0.01,1000],[Kd1, Kd1],'Color','red','LineStyle','--')
hold on
A=max((mean(ProbS(2:size(Kd2_list,2),:),2)/Tnum)/pA./(1-mean(ProbS(2:size(Kd2_list,2),:),2)/Tnum))
scatter(Kd2_list(2:size(Kd2_list,2)),(mean(ProbS(2:size(Kd2_list,2),:),2)/Tnum)/pA./(1-mean(ProbS(2:size(Kd2_list,2),:),2)/Tnum),'o','blue')
hold on
title(type+" / Kd_2 vs Kd_eff")
xlabel("log Kd2")
ylabel("log Kd_eff")
axis([0.01 1000 0 10000+0.1])
saveas(gcf,"Figure\"+Project_title+"___kd2vskdeff.png")

disp("save done")


%%
function ProbS_column = par_Metropolis_withW(Project_title,type,Tnum,Kd1,Kd2,Kd2_eff,pA,TestTime, MCMC_num)

ProbS_column=zeros(1,TestTime);

%disp("start simulation for kD2="+ string(kD2))
sys = Init_AT_System(type,Tnum);
for t=1:TestTime
    if rem(t,2^10)==0
        %disp("Done simulation for kD2="+ string(kD2_list(i))+" & t="+string(t/2^10)+" th 2^10")
    end
    sys = Reinitialize(sys);
    for j=1:MCMC_num
        sys = Metropolis_withW(sys,Kd1,Kd2,Kd2_eff,pA);
    end
    ProbS_column(1,t)=CalculateBindingNum(sys);
end
sys_model=sys;
save("Data\"+Project_title+"_\"+"Kd1_"+string(Kd1)+"___Kd2_"+string(Kd2)+".mat",'TestTime','Kd2','Kd1','Kd2_eff','pA','sys_model','Tnum','type','ProbS_column')

end

%%
function d=temp(Kd2, pA)
    d=1-roots([1 -2-Kd2/pA/8 1]);
    d=d(2);
end

%%
function ProbS_column = par_Metropolis_withW_temp(Project_title,type,Tnum,Kd1,Kd2,Kd2_eff,pA,TestTime, MCMC_num)

ProbS_column=zeros(1,TestTime);
d=temp(Kd2,pA);
pA=d*pA
%disp("start simulation for kD2="+ string(kD2))
sys = Init_AT_System(type,Tnum);
for t=1:TestTime
    if rem(t,2^10)==0
        %disp("Done simulation for kD2="+ string(kD2_list(i))+" & t="+string(t/2^10)+" th 2^10")
    end
    sys = Reinitialize(sys);
    for j=1:MCMC_num
        sys = Metropolis_withW(sys,Kd1,Kd2,Kd2_eff,pA);
    end
    ProbS_column(1,t)=CalculateBindingNum(sys);
end
sys_model=sys;
save("Data\"+Project_title+"_\"+"Kd1_"+string(Kd1)+"___Kd2_"+string(Kd2)+".mat",'TestTime','Kd2','Kd1','Kd2_eff','pA','sys_model','Tnum','type','ProbS_column')

end
