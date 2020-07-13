% Figure drawing v0.0 - 200701
%%
clc; clear;
cmap=hsv(10);
%% system type-wise plotting
data1name="triangularSphere2D"
datapath="Data\"+data1name+"_Tnum100_kD1_0.5_pA_1"+".mat"
load(datapath)
semilogx(kD2_list(2:size(kD2_list,2)),mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum,'-o','Color',cmap(1,:),'MarkerEdgeColor',cmap(1,:))
hold on


data2name="square2D"
datapath="Data\"+data2name+"_Tnum100_kD1_0.5_pA_1"+".mat"
load(datapath)
semilogx(kD2_list(2:size(kD2_list,2)),mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum,'-o','Color',cmap(2,:),'MarkerEdgeColor',cmap(2,:))
hold on


data3name="cubicSphere2D"
datapath="Data\"+data3name+"_Tnum100_kD1_0.5_pA_1"+".mat"
load(datapath)
semilogx(kD2_list(2:size(kD2_list,2)),mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum,'-o','Color',cmap(3,:),'MarkerEdgeColor',cmap(3,:))
hold on

data4name="linear1D"
datapath="Data\"+data4name+"_Tnum100_kD1_0.5_pA_1"+".mat"
load(datapath)
semilogx(kD2_list(2:size(kD2_list,2)),mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum,'-o','Color',cmap(4,:),'MarkerEdgeColor',cmap(4,:))
hold on

data5name="circular1D"
datapath="Data\"+data4name+"_Tnum100_kD1_0.5_pA_1"+".mat"
load(datapath)
semilogx(kD2_list(2:size(kD2_list,2)),mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum,'-o','Color',cmap(5,:),'MarkerEdgeColor',cmap(5,:))
hold on

semilogx([0.01,1000],[kD1, kD1],'Color','blue','LineStyle','--')
title("Kd_2 vs bounding %")
xlabel("log(Kd2)")
ylabel("bounding %")

legend([data1name data2name data3name data4name data5name "Control"])
a=string(datetime);
saveas(gcf,"Figure/"+a+".fig")
saveas(gcf,"Figure/"+a+".png")


%%
%%
data1name="triangularSphere2D"
datapath="Data\"+data1name+"_Tnum100_kD1_0.5_pA_1"+".mat"
load(datapath)
loglog(kD2_list(2:size(kD2_list,2)),(mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum)/pA./(1-mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum),'-o','Color',cmap(1,:),'MarkerEdgeColor',cmap(1,:))
hold on

data2name="square2D"
datapath="Data\"+data2name+"_Tnum100_kD1_0.5_pA_1"+".mat"
load(datapath)
loglog(kD2_list(2:size(kD2_list,2)),(mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum)/pA./(1-mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum),'-o','Color',cmap(2,:),'MarkerEdgeColor',cmap(2,:))
hold on

data3name="cubicSphere2D"
datapath="Data\"+data3name+"_Tnum100_kD1_0.5_pA_1"+".mat"
load(datapath)
loglog(kD2_list(2:size(kD2_list,2)),(mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum)/pA./(1-mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum),'-o','Color',cmap(3,:),'MarkerEdgeColor',cmap(3,:))
hold on

data4name="linear1D"
datapath="Data\"+data4name+"_Tnum100_kD1_0.5_pA_1"+".mat"
load(datapath)
loglog(kD2_list(2:size(kD2_list,2)),(mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum)/pA./(1-mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum),'-o','Color',cmap(4,:),'MarkerEdgeColor',cmap(4,:))
hold on

data5name="circular1D"
datapath="Data\"+data5name+"_Tnum100_kD1_0.5_pA_1"+".mat"
load(datapath)
loglog(kD2_list(2:size(kD2_list,2)),(mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum)/pA./(1-mean(ProbS(2:size(kD2_list,2),:),2)/sys.Tnum),'-o','Color',cmap(5,:),'MarkerEdgeColor',cmap(5,:))
hold on

semilogx([0.01,1000],[kD1, kD1],'Color','blue','LineStyle','--')
title("Kd_2 vs Kd_eff")
xlabel("log(Kd2)")
ylabel("log(Kd_eff)")

legend([data1name data2name data3name data4name data5name "Control"])
a=string(datetime);
saveas(gcf,"Figure/"+a+".fig")
saveas(gcf,"Figure/"+a+".png")

