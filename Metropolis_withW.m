function sys = Metropolis_withW(sys, kD1, kD2, kD2_eff, pA)
numIters = 2^10 * numel(sys.Tnum);

p_iter2=ceil(size(sys.W_relation,1)/sys.Tnum);
p1=min(pA/kD1,1);
p2=min(kD1/pA,1);
p3=min(1/kD2_eff,1);
p4=min(kD2_eff,1);
for iter1 = 1 : numIters
    % Pick a random target
    Ind = randi(sys.Tnum);
    % convert target
    if all(sys.W(sys.T2W{Ind})==0)
        if sys.T(Ind)==0
            if rand()<p1;
                sys.T(Ind)=1;
            end
        elseif sys.T(Ind)==1
            if rand()<p2;
                sys.T(Ind)=0;
            end
        end
    end
    
    for iter2 = 1 : p_iter2
        Ind = randi(size(sys.W_relation,1));
        if(sys.T(sys.W_relation(Ind,1)) & sys.T(sys.W_relation(Ind,2)))
            if sys.W(Ind)==0
                if rand()<p3;
                    sys.W(Ind)=1;
                end
            elseif sys.W(Ind)==1
                if rand()<p4;
                    sys.W(Ind)=0;
                end
            end
        end
    end
end