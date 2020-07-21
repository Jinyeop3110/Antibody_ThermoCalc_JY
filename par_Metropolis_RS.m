function ProbS_column = par_Metropolis_RS(Project_title,type,L,density,Kd1,Kd2,Kd2_eff,pA,TestTime, MCMC_num, WperT)


ProbS_column=zeros(1,TestTime);



%disp("start simulation for kD2="+ string(kD2))

sys = Init_RandomSurface_AT_System(type,L,density, WperT);

for t=1:TestTime

    if rem(t,2^10)==0

        %disp("Done simulation for kD2="+ string(Kd2_list(i))+" & t="+string(t/2^10)+" th 2^10")

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