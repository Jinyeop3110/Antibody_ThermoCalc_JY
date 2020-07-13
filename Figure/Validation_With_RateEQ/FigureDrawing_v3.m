% Figure drawing v3 - 200712
%% For validation with theory
clc; clear;
cmap=hsv(10);
%% 
j=0;
lgd=[];

j=j+1;
Tnum=2;
datapath = "Data\0711_linear1D_Tnum"+int2str(Tnum)+"_Nchange.mat";
dataname="Mycode Tnum="+int2str(Tnum);
lgd=[lgd dataname];
load(datapath)
semilogx(kD2_list(2:size(kD2_list,2)),mean(ProbS(2:size(kD2_list,2),:),2)/sys_model.Tnum,'-o','Color',cmap(j,:),'MarkerEdgeColor',cmap(j,:))
hold on

j=j+5;
dataname="Theory Tnum="+int2str(Tnum);
lgd=[lgd dataname];
t=pA/kD1;
s=1./kD2_eff_list;
bp=(t+t^2+s*t^2)./(1+2*t+t^2+s*t^2);
semilogx(kD2_list(2:size(kD2_list,2)),bp(2:size(kD2_list,2)),'-','LineWidth',1.2,'Color',cmap(j,:),'MarkerEdgeColor',cmap(j,:))
hold on;

semilogx([0.5*min(kD2_list(2:size(kD2_list,2))),2*max(kD2_list(2:size(kD2_list,2)))],[pA/kD1, pA/kD1],'LineWidth',1.2,'Color','blue','LineStyle','--')
title("Kd_2 vs bounding #")
xlabel("log(Kd2)")
ylabel("bounding #")
lgd=[lgd "Control"];

legend(lgd)
a=[sprintf('%s', datestr(now,'mm-dd-yyyy HH-MM_')) int2str(randi(500))];
saveas(gcf,"Figure\fig_"+a+".fig")
saveas(gcf,"Figure\fig_"+a+".png")

%%
figure()
j=0;
lgd=[];

j=j+1;
Tnum=3;
datapath = "Data\0711_linear1D_Tnum"+int2str(Tnum)+"_Nchange.mat";
dataname="Mycode Tnum="+int2str(Tnum);
lgd=[lgd dataname];
load(datapath)
semilogx(kD2_list(2:size(kD2_list,2)),mean(ProbS(2:size(kD2_list,2),:),2)/sys_model.Tnum,'-o','LineWidth',1.2,'Color',cmap(j,:),'MarkerEdgeColor',cmap(j,:))
hold on

j=j+5;
dataname="Theory Tnum="+int2str(Tnum);
lgd=[lgd dataname];
t=pA/kD1;
s=1./kD2_eff_list;
bp=(3*t+6*t^2+4*t^2*s+3*t^3+6*t^3*s+3*t^3*s.^2)./(1+3*t+3*t^2+2*t^2*s+t^3+2*t^3*s+t^3*s.^2)/3;
semilogx(kD2_list(2:size(kD2_list,2)),bp(2:size(kD2_list,2)),'-','Color',cmap(j,:),'LineWidth',1.2,'MarkerEdgeColor',cmap(j,:))
hold on;

semilogx([0.5*min(kD2_list(2:size(kD2_list,2))),2*max(kD2_list(2:size(kD2_list,2)))],[pA/kD1, pA/kD1],'Color','blue','LineStyle','--')
title("Kd_2 vs bounding #")
xlabel("log(Kd2)")
ylabel("bounding #")
lgd=[lgd "Control"];

legend(lgd)
a=[sprintf('%s', datestr(now,'mm-dd-yyyy HH-MM_')) int2str(randi(500))];
saveas(gcf,"Figure\fig_"+a+".fig")
saveas(gcf,"Figure\fig_"+a+".png")

%%
%% 
lgd=[];
j=0;

for j=1:8
    j=j+1;
Tnum=Tnum_list(j);
datapath = "Data\0711_square2D_Tnum"+int2str(Tnum)+"_Nchange.mat";
dataname="Tnum="+int2str(Tnum);
lgd=[lgd dataname];
load(datapath)
loglog(kD2_list(2:size(kD2_list,2)),pA*(1-mean(ProbS(2:size(kD2_list,2),:),2)/sys_model.Tnum)./(mean(ProbS(2:size(kD2_list,2),:),2)/sys_model.Tnum),'-o','Color',cmap(j,:),'MarkerEdgeColor',cmap(2,:))
hold on
end



semilogx([0.5*min(kD2_list(2:size(kD2_list,2))),2*max(kD2_list(2:size(kD2_list,2)))],[kD1, kD1],'Color','blue','LineStyle','--')
title("Kd_2 vs Kd_eff")
xlabel("log(Kd2)")
ylabel("log(Kd_eff)")
legend(lgd)

legend(lgd)
a=[sprintf('%s', datestr(now,'mm-dd-yyyy HH-MM_')) int2str(randi(500))];
saveas(gcf,"Figure\fig_"+a+".fig")
saveas(gcf,"Figure\fig_"+a+".png")
