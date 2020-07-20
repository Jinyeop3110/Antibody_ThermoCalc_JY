function sys = Init_RandomSurface_AT_System(type, L, density)

% RandomSurface AT system initiation , Jul-20-2020

%   Copyright 2020 Jinyeop Song.



sys=AT_System;



%% randomUniformFlat2D
if(type=="randomUniformFlat2D")
    p=0.5;
    sys.type="randomUniformFlat2D";
    sys.Tnum=floor((L)*density);
    sys.T=(rand(1,sys.Tnum)>p);
    sys.T_position=struct;
    
    sys.T_position.x=(L^0.5)*rand(sys.Tnum,1);
    sys.T_position.y=(L^0.5)*rand(sys.Tnum,1);
    sys.T_position.z=0;
    
    sys.W_relation=[];
    sys.T2W=cell(1,sys.Tnum);
    wi=0;
    
    for 
    
    for i=1:sys.Tnum
        for j=(i+1):sys.Tnum
            if(DistanceBtwTwoTargets(sys,i,j)<1.01)
                wi=wi+1;
                sys.W_relation = [sys.W_relation ; [i,j]];
                sys.T2W{i}=[sys.T2W{i}, wi];
                sys.T2W{j}=[sys.T2W{j}, wi];
            end
        end
    end
    sys.W=zeros(1,size(sys.W_relation,1));
end

%% randomUniformSphere2D
if(type=="randomUniformSphere2D")
    p=0.5;
    sys.type="randomUniformSphere2D";
    sys.Tnum=floor((L)*density);
    sys.T=(rand(1,sys.Tnum)>p);
    sys.T_position=struct;
    
    
    theta=arccos(1-2*rand(sys.Tnum,1))
    phi=2*3.141592*rand(sys.Tnum,1)
    radius=(L/4/3.141592)^0.5
    sys.T_position.x=radius*sin(theta)*cos(phi);
    sys.T_position.y=radius*sin(theta)*sin(phi);
    sys.T_position.z=cos(theta);
    
    sys.W_relation=[];
    sys.T2W=cell(1,sys.Tnum);
    wi=0;
    
    for 
    
    for i=1:sys.Tnum
        for j=(i+1):sys.Tnum
            if(DistanceBtwTwoTargets(sys,i,j)<1.01)
                wi=wi+1;
                sys.W_relation = [sys.W_relation ; [i,j]];
                sys.T2W{i}=[sys.T2W{i}, wi];
                sys.T2W{j}=[sys.T2W{j}, wi];
            end
        end
    end
    sys.W=zeros(1,size(sys.W_relation,1));
end
%% etc

%% randomQuasiFlat2D
if(type=="randomQuasiFlat2D")
    p=0.5;
    sys.type="randomQuasiFlat2D";
    sys.Tnum=floor((L)*density);
    sys.T=(rand(1,sys.Tnum)>p);
    sys.T_position=struct;
    
    p = haltonset(2,'Skip',1e3,'Leap',1e2)
    p = scramble(p,'RR2')
    X0 = net(p,sys.Tnum);

    sys.T_position.x=(L^0.5)*X0(:,1);
    sys.T_position.y=(L^0.5)*X0(:,2);
    sys.T_position.z=0;
    
    sys.W_relation=[];
    sys.T2W=cell(1,sys.Tnum);
    wi=0;
    
    for 
    
    for i=1:sys.Tnum
        for j=(i+1):sys.Tnum
            if(DistanceBtwTwoTargets(sys,i,j)<1.01)
                wi=wi+1;
                sys.W_relation = [sys.W_relation ; [i,j]];
                sys.T2W{i}=[sys.T2W{i}, wi];
                sys.T2W{j}=[sys.T2W{j}, wi];
            end
        end
    end
    sys.W=zeros(1,size(sys.W_relation,1));
end

%% randomQuasiSphere2D
if(type=="randomQuasiSphere2D")
    p=0.5;
    sys.type="randomQuasiSphere2D";
    sys.Tnum=floor((L)*density);
    sys.T=(rand(1,sys.Tnum)>p);
    sys.T_position=struct;
    
    
    p = haltonset(2,'Skip',1e3,'Leap',1e2)
    p = scramble(p,'RR2')
    X0 = net(p,sys.Tnum);

    sys.T_position.x=(L^0.5)*X0(:,1);
    
    theta=arccos(1-2*X0(:,1))
    phi=2*3.141592*X0(:,1)
    radius=(L/4/3.141592)^0.5
    sys.T_position.x=radius*sin(theta)*cos(phi);
    sys.T_position.y=radius*sin(theta)*sin(phi);
    sys.T_position.z=cos(theta);
    
    sys.W_relation=[];
    sys.T2W=cell(1,sys.Tnum);
    wi=0;
    
    for 
    
    for i=1:sys.Tnum
        for j=(i+1):sys.Tnum
            if(DistanceBtwTwoTargets(sys,i,j)<1.01)
                wi=wi+1;
                sys.W_relation = [sys.W_relation ; [i,j]];
                sys.T2W{i}=[sys.T2W{i}, wi];
                sys.T2W{j}=[sys.T2W{j}, wi];
            end
        end
    end
    sys.W=zeros(1,size(sys.W_relation,1));
end
%% etc

if(type=='etc')

end



disp(sys)

end
