% Figure drawing v0.0 - 200706
%%
clc; clear;
cmap=hsv(10);
%% 
data1name="cubicSphere2D"
datapath="Data\0705_cubicSphere2D_Tnum98.mat"
load(datapath)
semilogx(kD2_list(2:size(kD2_list,2)),mean(ProbS(2:size(kD2_list,2),:),2)/sys_model.Tnum,'-o','Color',cmap(1,:),'MarkerEdgeColor',cmap(1,:))
hold on

data2name="triangularSphere2D"
datapath="Data\0705_triangularSphere2D_Tnum98.mat"
load(datapath)
semilogx(kD2_list(2:size(kD2_list,2)),mean(ProbS(2:size(kD2_list,2),:),2)/sys_model.Tnum,'-o','Color',cmap(2,:),'MarkerEdgeColor',cmap(2,:))
hold on

semilogx([0.5*min(kD2_list(2:size(kD2_list,2))),2*max(kD2_list(2:size(kD2_list,2)))],[pA/kD1, pA/kD1],'Color','blue','LineStyle','--')
title("Kd_2 vs Kd_eff")
xlabel("log(Kd2)")
ylabel("bounding #")

legend([data1name data2name "Control"])

%%
%%
figure()
data1name="cubicSphere2D"
datapath="Data\0705_cubicSphere2D_Tnum98.mat"
load(datapath)
loglog(kD2_list(2:size(kD2_list,2)),pA*(1-mean(ProbS(2:size(kD2_list,2),:),2)/sys_model.Tnum)./(mean(ProbS(2:size(kD2_list,2),:),2)/sys_model.Tnum),'-o','Color',cmap(1,:),'MarkerEdgeColor',cmap(1,:))
hold on

data2name="triangularSphere2D"
datapath="Data\0705_triangularSphere2D_Tnum98.mat"
load(datapath)
loglog(kD2_list(2:size(kD2_list,2)),pA*(1-mean(ProbS(2:size(kD2_list,2),:),2)/sys_model.Tnum)./(mean(ProbS(2:size(kD2_list,2),:),2)/sys_model.Tnum),'-o','Color',cmap(2,:),'MarkerEdgeColor',cmap(2,:))
hold on

semilogx([0.5*min(kD2_list(2:size(kD2_list,2))),2*max(kD2_list(2:size(kD2_list,2)))],[kD1, kD1],'Color','blue','LineStyle','--')
title("Kd_2 vs Kd_eff")
xlabel("log(Kd2)")
ylabel("log(Kd_eff)")

legend([data1name data2name "Control"])
