clear all;
close all;
%------------- PSO code
    %% Parameter
    nSW=50;  % number of particles
    nPam=4; % dimension; number of design variable
    nIter_max=100;  % maximum of iterations
    min_fit_criteria=0; % acceptable minimum value of fitness
    inertia_weight=0.72;
    acc_coef1=2;
    acc_coef2=2;
    ub=600;
    lb=-600;
    ubt=600;
    lbt=-600;
    %  -----------------------------------------------------
    %           Initializing
        for tb=1:1000

    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    %---------------------------------------------------
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        Vel=inertia_weight*Vel+acc_coef1*rand(nSW,nPam).*(Pbest-X)+acc_coef2*rand(nSW,nPam).*(ones(nSW,1)*Gbest-X);
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    psoGriewangk4(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
        end
fprintf(' Ýter=%d ObjVal=%g\n',var(psoGriewangk4),mean(psoGriewangk4))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('psoGriewangk4','psoGriewangk4')





clear all;
close all;
%------------- PSO code
%% Parameter
nSW=50;  % number of particles
nPam=4; % dimension; number of design variable
nIter_max=100;  % maximum of iterations
min_fit_criteria=0; % acceptable minimum value of fitness
inertia_weight=0.72;
acc_coef1=eye(nSW);
acc_coef2=eye(nSW);
ub=600;
lb=-600;
ubt=600;
lbt=-600;

%---------------------------------------------------
%           fuzzy Parameter
nbn=newfis('fpso');
nbn=addvar(nbn,'input','iter',[0 nIter_max]);
nbn=addmf(nbn,'input',1,'start','gaussmf',[.25*nIter_max 0]);
nbn=addmf(nbn,'input',1,'mid2','gaussmf',[.25*nIter_max .5*nIter_max]);
nbn=addmf(nbn,'input',1,'end','gaussmf',[.25*nIter_max nIter_max]);

nbn=addvar(nbn,'input','cfval',[0 1]);
nbn=addmf(nbn,'input',2,'bad','gaussmf',[.25 0]);
nbn=addmf(nbn,'input',2,'normal','gaussmf',[.25 .5]);
nbn=addmf(nbn,'input',2,'excellent','gaussmf',[.25 1]);

nbn=addvar(nbn,'output','c1',[0 3]);
nbn=addmf(nbn,'output',1,'lowc1','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',1,'medium1c1','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',1,'medium2c1','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',1,'medium3c1','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',1,'highc1','trimf',[2 2.5 3]);

nbn=addvar(nbn,'output','c2',[0 3]);
nbn=addmf(nbn,'output',2,'lowc2','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',2,'medium1c2','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',2,'medium2c2','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',2,'medium3c2','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',2,'highc2','trimf',[2 2.5 3]);

ruleList=[ ...
    1 1 3 3 1 2
    1 2 2 4 1 2
    1 3 1 5 1 2
    2 1 4 2 1 2
    2 2 3 3 1 2
    2 3 2 4 1 2
    3 1 5 1 1 2
    3 2 4 2 1 2
    3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 4 2 1 2
%     1 3 5 1 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 1 5 1 2
%     3 2 2 4 1 2
%     3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 1 5 1 2
%     1 2 2 4 1 2
%     1 3 3 3 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 3 3 1 2
%     3 2 4 2 1 2
%     3 3 5 1 1 2];
% ruleList=[ ...
%     1 1 5 1 1 2
%     1 2 4 2 1 2
%     1 3 3 3 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 3 3 1 2
%     3 2 2 4 1 2
%     3 3 1 5 1 2];

nbn=addrule(nbn,ruleList);

%  -----------------------------------------------------
%           Initializing
for tb=1:1000
    
    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        mincf=min(F);
        maxcf=max(F);
        normf=(F-mincf)/(maxcf-mincf);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        ansc=evalfis([i*ones(nSW,1) normf],nbn);
        acc_coef1(acc_coef1~=0)=ansc(:,1);
        acc_coef2(acc_coef2~=0)=ansc(:,2);
        Vel=inertia_weight*Vel+acc_coef1*(rand(nSW,nPam).*(Pbest-X))+acc_coef2*(rand(nSW,nPam).*(ones(nSW,1)*Gbest-X));
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    fuzpsoGriewangk4r1(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
end
fprintf(' Ýter=%d ObjVal=%g\n',var(fuzpsoGriewangk4r1),mean(fuzpsoGriewangk4r1))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('fuzpsoGriewangk4r1','fuzpsoGriewangk4r1')






clear all;
close all;
%------------- PSO code
%% Parameter
nSW=50;  % number of particles
nPam=4; % dimension; number of design variable
nIter_max=100;  % maximum of iterations
min_fit_criteria=0; % acceptable minimum value of fitness
inertia_weight=0.72;
acc_coef1=eye(nSW);
acc_coef2=eye(nSW);
ub=600;
lb=-600;
ubt=600;
lbt=-600;

%---------------------------------------------------
%           fuzzy Parameter
nbn=newfis('fpso');
nbn=addvar(nbn,'input','iter',[0 nIter_max]);
nbn=addmf(nbn,'input',1,'start','gaussmf',[.25*nIter_max 0]);
nbn=addmf(nbn,'input',1,'mid2','gaussmf',[.25*nIter_max .5*nIter_max]);
nbn=addmf(nbn,'input',1,'end','gaussmf',[.25*nIter_max nIter_max]);

nbn=addvar(nbn,'input','cfval',[0 1]);
nbn=addmf(nbn,'input',2,'bad','gaussmf',[.25 0]);
nbn=addmf(nbn,'input',2,'normal','gaussmf',[.25 .5]);
nbn=addmf(nbn,'input',2,'excellent','gaussmf',[.25 1]);

nbn=addvar(nbn,'output','c1',[0 3]);
nbn=addmf(nbn,'output',1,'lowc1','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',1,'medium1c1','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',1,'medium2c1','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',1,'medium3c1','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',1,'highc1','trimf',[2 2.5 3]);

nbn=addvar(nbn,'output','c2',[0 3]);
nbn=addmf(nbn,'output',2,'lowc2','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',2,'medium1c2','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',2,'medium2c2','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',2,'medium3c2','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',2,'highc2','trimf',[2 2.5 3]);

% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 2 4 1 2
%     1 3 1 5 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 5 1 1 2
%     3 2 4 2 1 2
%     3 3 3 3 1 2];
ruleList=[ ...
    1 1 3 3 1 2
    1 2 4 2 1 2
    1 3 5 1 1 2
    2 1 2 4 1 2
    2 2 3 3 1 2
    2 3 4 2 1 2
    3 1 1 5 1 2
    3 2 2 4 1 2
    3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 1 5 1 2
%     1 2 2 4 1 2
%     1 3 3 3 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 3 3 1 2
%     3 2 4 2 1 2
%     3 3 5 1 1 2];
% ruleList=[ ...
%     1 1 5 1 1 2
%     1 2 4 2 1 2
%     1 3 3 3 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 3 3 1 2
%     3 2 2 4 1 2
%     3 3 1 5 1 2];

nbn=addrule(nbn,ruleList);

%  -----------------------------------------------------
%           Initializing
for tb=1:1000
    
    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        mincf=min(F);
        maxcf=max(F);
        normf=(F-mincf)/(maxcf-mincf);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        ansc=evalfis([i*ones(nSW,1) normf],nbn);
        acc_coef1(acc_coef1~=0)=ansc(:,1);
        acc_coef2(acc_coef2~=0)=ansc(:,2);
        Vel=inertia_weight*Vel+acc_coef1*(rand(nSW,nPam).*(Pbest-X))+acc_coef2*(rand(nSW,nPam).*(ones(nSW,1)*Gbest-X));
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    fuzpsoGriewangk4r2(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
end
fprintf(' Ýter=%d ObjVal=%g\n',var(fuzpsoGriewangk4r2),mean(fuzpsoGriewangk4r2))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('fuzpsoGriewangk4r2','fuzpsoGriewangk4r2')






clear all;
close all;
%------------- PSO code
%% Parameter
nSW=50;  % number of particles
nPam=4; % dimension; number of design variable
nIter_max=100;  % maximum of iterations
min_fit_criteria=0; % acceptable minimum value of fitness
inertia_weight=0.72;
acc_coef1=eye(nSW);
acc_coef2=eye(nSW);
ub=600;
lb=-600;
ubt=600;
lbt=-600;

%---------------------------------------------------
%           fuzzy Parameter
nbn=newfis('fpso');
nbn=addvar(nbn,'input','iter',[0 nIter_max]);
nbn=addmf(nbn,'input',1,'start','gaussmf',[.25*nIter_max 0]);
nbn=addmf(nbn,'input',1,'mid2','gaussmf',[.25*nIter_max .5*nIter_max]);
nbn=addmf(nbn,'input',1,'end','gaussmf',[.25*nIter_max nIter_max]);

nbn=addvar(nbn,'input','cfval',[0 1]);
nbn=addmf(nbn,'input',2,'bad','gaussmf',[.25 0]);
nbn=addmf(nbn,'input',2,'normal','gaussmf',[.25 .5]);
nbn=addmf(nbn,'input',2,'excellent','gaussmf',[.25 1]);

nbn=addvar(nbn,'output','c1',[0 3]);
nbn=addmf(nbn,'output',1,'lowc1','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',1,'medium1c1','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',1,'medium2c1','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',1,'medium3c1','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',1,'highc1','trimf',[2 2.5 3]);

nbn=addvar(nbn,'output','c2',[0 3]);
nbn=addmf(nbn,'output',2,'lowc2','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',2,'medium1c2','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',2,'medium2c2','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',2,'medium3c2','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',2,'highc2','trimf',[2 2.5 3]);

% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 2 4 1 2
%     1 3 1 5 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 5 1 1 2
%     3 2 4 2 1 2
%     3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 4 2 1 2
%     1 3 5 1 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 1 5 1 2
%     3 2 2 4 1 2
%     3 3 3 3 1 2];
ruleList=[ ...
    1 1 1 5 1 2
    1 2 2 4 1 2
    1 3 3 3 1 2
    2 1 2 4 1 2
    2 2 3 3 1 2
    2 3 4 2 1 2
    3 1 3 3 1 2
    3 2 4 2 1 2
    3 3 5 1 1 2];
% ruleList=[ ...
%     1 1 5 1 1 2
%     1 2 4 2 1 2
%     1 3 3 3 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 3 3 1 2
%     3 2 2 4 1 2
%     3 3 1 5 1 2];

nbn=addrule(nbn,ruleList);

%  -----------------------------------------------------
%           Initializing
for tb=1:1000
    
    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        mincf=min(F);
        maxcf=max(F);
        normf=(F-mincf)/(maxcf-mincf);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        ansc=evalfis([i*ones(nSW,1) normf],nbn);
        acc_coef1(acc_coef1~=0)=ansc(:,1);
        acc_coef2(acc_coef2~=0)=ansc(:,2);
        Vel=inertia_weight*Vel+acc_coef1*(rand(nSW,nPam).*(Pbest-X))+acc_coef2*(rand(nSW,nPam).*(ones(nSW,1)*Gbest-X));
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    fuzpsoGriewangk4r3(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
end
fprintf(' Ýter=%d ObjVal=%g\n',var(fuzpsoGriewangk4r3),mean(fuzpsoGriewangk4r3))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('fuzpsoGriewangk4r3','fuzpsoGriewangk4r3')





clear all;
close all;
%------------- PSO code
%% Parameter
nSW=50;  % number of particles
nPam=4; % dimension; number of design variable
nIter_max=100;  % maximum of iterations
min_fit_criteria=0; % acceptable minimum value of fitness
inertia_weight=0.72;
acc_coef1=eye(nSW);
acc_coef2=eye(nSW);
ub=600;
lb=-600;
ubt=600;
lbt=-600;

%---------------------------------------------------
%           fuzzy Parameter
nbn=newfis('fpso');
nbn=addvar(nbn,'input','iter',[0 nIter_max]);
nbn=addmf(nbn,'input',1,'start','gaussmf',[.25*nIter_max 0]);
nbn=addmf(nbn,'input',1,'mid2','gaussmf',[.25*nIter_max .5*nIter_max]);
nbn=addmf(nbn,'input',1,'end','gaussmf',[.25*nIter_max nIter_max]);

nbn=addvar(nbn,'input','cfval',[0 1]);
nbn=addmf(nbn,'input',2,'bad','gaussmf',[.25 0]);
nbn=addmf(nbn,'input',2,'normal','gaussmf',[.25 .5]);
nbn=addmf(nbn,'input',2,'excellent','gaussmf',[.25 1]);

nbn=addvar(nbn,'output','c1',[0 3]);
nbn=addmf(nbn,'output',1,'lowc1','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',1,'medium1c1','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',1,'medium2c1','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',1,'medium3c1','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',1,'highc1','trimf',[2 2.5 3]);

nbn=addvar(nbn,'output','c2',[0 3]);
nbn=addmf(nbn,'output',2,'lowc2','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',2,'medium1c2','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',2,'medium2c2','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',2,'medium3c2','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',2,'highc2','trimf',[2 2.5 3]);

% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 2 4 1 2
%     1 3 1 5 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 5 1 1 2
%     3 2 4 2 1 2
%     3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 4 2 1 2
%     1 3 5 1 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 1 5 1 2
%     3 2 2 4 1 2
%     3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 1 5 1 2
%     1 2 2 4 1 2
%     1 3 3 3 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 3 3 1 2
%     3 2 4 2 1 2
%     3 3 5 1 1 2];
ruleList=[ ...
    1 1 5 1 1 2
    1 2 4 2 1 2
    1 3 3 3 1 2
    2 1 4 2 1 2
    2 2 3 3 1 2
    2 3 2 4 1 2
    3 1 3 3 1 2
    3 2 2 4 1 2
    3 3 1 5 1 2];

nbn=addrule(nbn,ruleList);

%  -----------------------------------------------------
%           Initializing
for tb=1:1000
    
    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        mincf=min(F);
        maxcf=max(F);
        normf=(F-mincf)/(maxcf-mincf);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        ansc=evalfis([i*ones(nSW,1) normf],nbn);
        acc_coef1(acc_coef1~=0)=ansc(:,1);
        acc_coef2(acc_coef2~=0)=ansc(:,2);
        Vel=inertia_weight*Vel+acc_coef1*(rand(nSW,nPam).*(Pbest-X))+acc_coef2*(rand(nSW,nPam).*(ones(nSW,1)*Gbest-X));
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    fuzpsoGriewangk4r4(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
end
fprintf(' Ýter=%d ObjVal=%g\n',var(fuzpsoGriewangk4r4),mean(fuzpsoGriewangk4r4))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('fuzpsoGriewangk4r4','fuzpsoGriewangk4r4')





clear all;
close all;
%------------- PSO code
    %% Parameter
    nSW=50;  % number of particles
    nPam=8; % dimension; number of design variable
    nIter_max=100;  % maximum of iterations
    min_fit_criteria=0; % acceptable minimum value of fitness
    inertia_weight=0.72;
    acc_coef1=2;
    acc_coef2=2;
    ub=600;
    lb=-600;
    ubt=600;
    lbt=-600;
    %  -----------------------------------------------------
    %           Initializing
        for tb=1:1000

    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    %---------------------------------------------------
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        Vel=inertia_weight*Vel+acc_coef1*rand(nSW,nPam).*(Pbest-X)+acc_coef2*rand(nSW,nPam).*(ones(nSW,1)*Gbest-X);
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    psoGriewangk8(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
        end
fprintf(' Ýter=%d ObjVal=%g\n',var(psoGriewangk8),mean(psoGriewangk8))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('psoGriewangk8','psoGriewangk8')





clear all;
close all;
%------------- PSO code
%% Parameter
nSW=50;  % number of particles
nPam=8; % dimension; number of design variable
nIter_max=100;  % maximum of iterations
min_fit_criteria=0; % acceptable minimum value of fitness
inertia_weight=0.72;
acc_coef1=eye(nSW);
acc_coef2=eye(nSW);
ub=600;
lb=-600;
ubt=600;
lbt=-600;

%---------------------------------------------------
%           fuzzy Parameter
nbn=newfis('fpso');
nbn=addvar(nbn,'input','iter',[0 nIter_max]);
nbn=addmf(nbn,'input',1,'start','gaussmf',[.25*nIter_max 0]);
nbn=addmf(nbn,'input',1,'mid2','gaussmf',[.25*nIter_max .5*nIter_max]);
nbn=addmf(nbn,'input',1,'end','gaussmf',[.25*nIter_max nIter_max]);

nbn=addvar(nbn,'input','cfval',[0 1]);
nbn=addmf(nbn,'input',2,'bad','gaussmf',[.25 0]);
nbn=addmf(nbn,'input',2,'normal','gaussmf',[.25 .5]);
nbn=addmf(nbn,'input',2,'excellent','gaussmf',[.25 1]);

nbn=addvar(nbn,'output','c1',[0 3]);
nbn=addmf(nbn,'output',1,'lowc1','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',1,'medium1c1','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',1,'medium2c1','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',1,'medium3c1','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',1,'highc1','trimf',[2 2.5 3]);

nbn=addvar(nbn,'output','c2',[0 3]);
nbn=addmf(nbn,'output',2,'lowc2','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',2,'medium1c2','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',2,'medium2c2','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',2,'medium3c2','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',2,'highc2','trimf',[2 2.5 3]);

ruleList=[ ...
    1 1 3 3 1 2
    1 2 2 4 1 2
    1 3 1 5 1 2
    2 1 4 2 1 2
    2 2 3 3 1 2
    2 3 2 4 1 2
    3 1 5 1 1 2
    3 2 4 2 1 2
    3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 4 2 1 2
%     1 3 5 1 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 1 5 1 2
%     3 2 2 4 1 2
%     3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 1 5 1 2
%     1 2 2 4 1 2
%     1 3 3 3 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 3 3 1 2
%     3 2 4 2 1 2
%     3 3 5 1 1 2];
% ruleList=[ ...
%     1 1 5 1 1 2
%     1 2 4 2 1 2
%     1 3 3 3 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 3 3 1 2
%     3 2 2 4 1 2
%     3 3 1 5 1 2];

nbn=addrule(nbn,ruleList);

%  -----------------------------------------------------
%           Initializing
for tb=1:1000
    
    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        mincf=min(F);
        maxcf=max(F);
        normf=(F-mincf)/(maxcf-mincf);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        ansc=evalfis([i*ones(nSW,1) normf],nbn);
        acc_coef1(acc_coef1~=0)=ansc(:,1);
        acc_coef2(acc_coef2~=0)=ansc(:,2);
        Vel=inertia_weight*Vel+acc_coef1*(rand(nSW,nPam).*(Pbest-X))+acc_coef2*(rand(nSW,nPam).*(ones(nSW,1)*Gbest-X));
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    fuzpsoGriewangk8r1(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
end
fprintf(' Ýter=%d ObjVal=%g\n',var(fuzpsoGriewangk8r1),mean(fuzpsoGriewangk8r1))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('fuzpsoGriewangk8r1','fuzpsoGriewangk8r1')






clear all;
close all;
%------------- PSO code
%% Parameter
nSW=50;  % number of particles
nPam=8; % dimension; number of design variable
nIter_max=100;  % maximum of iterations
min_fit_criteria=0; % acceptable minimum value of fitness
inertia_weight=0.72;
acc_coef1=eye(nSW);
acc_coef2=eye(nSW);
ub=600;
lb=-600;
ubt=600;
lbt=-600;

%---------------------------------------------------
%           fuzzy Parameter
nbn=newfis('fpso');
nbn=addvar(nbn,'input','iter',[0 nIter_max]);
nbn=addmf(nbn,'input',1,'start','gaussmf',[.25*nIter_max 0]);
nbn=addmf(nbn,'input',1,'mid2','gaussmf',[.25*nIter_max .5*nIter_max]);
nbn=addmf(nbn,'input',1,'end','gaussmf',[.25*nIter_max nIter_max]);

nbn=addvar(nbn,'input','cfval',[0 1]);
nbn=addmf(nbn,'input',2,'bad','gaussmf',[.25 0]);
nbn=addmf(nbn,'input',2,'normal','gaussmf',[.25 .5]);
nbn=addmf(nbn,'input',2,'excellent','gaussmf',[.25 1]);

nbn=addvar(nbn,'output','c1',[0 3]);
nbn=addmf(nbn,'output',1,'lowc1','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',1,'medium1c1','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',1,'medium2c1','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',1,'medium3c1','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',1,'highc1','trimf',[2 2.5 3]);

nbn=addvar(nbn,'output','c2',[0 3]);
nbn=addmf(nbn,'output',2,'lowc2','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',2,'medium1c2','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',2,'medium2c2','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',2,'medium3c2','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',2,'highc2','trimf',[2 2.5 3]);

% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 2 4 1 2
%     1 3 1 5 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 5 1 1 2
%     3 2 4 2 1 2
%     3 3 3 3 1 2];
ruleList=[ ...
    1 1 3 3 1 2
    1 2 4 2 1 2
    1 3 5 1 1 2
    2 1 2 4 1 2
    2 2 3 3 1 2
    2 3 4 2 1 2
    3 1 1 5 1 2
    3 2 2 4 1 2
    3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 1 5 1 2
%     1 2 2 4 1 2
%     1 3 3 3 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 3 3 1 2
%     3 2 4 2 1 2
%     3 3 5 1 1 2];
% ruleList=[ ...
%     1 1 5 1 1 2
%     1 2 4 2 1 2
%     1 3 3 3 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 3 3 1 2
%     3 2 2 4 1 2
%     3 3 1 5 1 2];

nbn=addrule(nbn,ruleList);

%  -----------------------------------------------------
%           Initializing
for tb=1:1000
    
    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        mincf=min(F);
        maxcf=max(F);
        normf=(F-mincf)/(maxcf-mincf);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        ansc=evalfis([i*ones(nSW,1) normf],nbn);
        acc_coef1(acc_coef1~=0)=ansc(:,1);
        acc_coef2(acc_coef2~=0)=ansc(:,2);
        Vel=inertia_weight*Vel+acc_coef1*(rand(nSW,nPam).*(Pbest-X))+acc_coef2*(rand(nSW,nPam).*(ones(nSW,1)*Gbest-X));
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    fuzpsoGriewangk8r2(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
end
fprintf(' Ýter=%d ObjVal=%g\n',var(fuzpsoGriewangk8r2),mean(fuzpsoGriewangk8r2))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('fuzpsoGriewangk8r2','fuzpsoGriewangk8r2')






clear all;
close all;
%------------- PSO code
%% Parameter
nSW=50;  % number of particles
nPam=8; % dimension; number of design variable
nIter_max=100;  % maximum of iterations
min_fit_criteria=0; % acceptable minimum value of fitness
inertia_weight=0.72;
acc_coef1=eye(nSW);
acc_coef2=eye(nSW);
ub=600;
lb=-600;
ubt=600;
lbt=-600;

%---------------------------------------------------
%           fuzzy Parameter
nbn=newfis('fpso');
nbn=addvar(nbn,'input','iter',[0 nIter_max]);
nbn=addmf(nbn,'input',1,'start','gaussmf',[.25*nIter_max 0]);
nbn=addmf(nbn,'input',1,'mid2','gaussmf',[.25*nIter_max .5*nIter_max]);
nbn=addmf(nbn,'input',1,'end','gaussmf',[.25*nIter_max nIter_max]);

nbn=addvar(nbn,'input','cfval',[0 1]);
nbn=addmf(nbn,'input',2,'bad','gaussmf',[.25 0]);
nbn=addmf(nbn,'input',2,'normal','gaussmf',[.25 .5]);
nbn=addmf(nbn,'input',2,'excellent','gaussmf',[.25 1]);

nbn=addvar(nbn,'output','c1',[0 3]);
nbn=addmf(nbn,'output',1,'lowc1','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',1,'medium1c1','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',1,'medium2c1','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',1,'medium3c1','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',1,'highc1','trimf',[2 2.5 3]);

nbn=addvar(nbn,'output','c2',[0 3]);
nbn=addmf(nbn,'output',2,'lowc2','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',2,'medium1c2','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',2,'medium2c2','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',2,'medium3c2','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',2,'highc2','trimf',[2 2.5 3]);

% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 2 4 1 2
%     1 3 1 5 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 5 1 1 2
%     3 2 4 2 1 2
%     3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 4 2 1 2
%     1 3 5 1 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 1 5 1 2
%     3 2 2 4 1 2
%     3 3 3 3 1 2];
ruleList=[ ...
    1 1 1 5 1 2
    1 2 2 4 1 2
    1 3 3 3 1 2
    2 1 2 4 1 2
    2 2 3 3 1 2
    2 3 4 2 1 2
    3 1 3 3 1 2
    3 2 4 2 1 2
    3 3 5 1 1 2];
% ruleList=[ ...
%     1 1 5 1 1 2
%     1 2 4 2 1 2
%     1 3 3 3 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 3 3 1 2
%     3 2 2 4 1 2
%     3 3 1 5 1 2];

nbn=addrule(nbn,ruleList);

%  -----------------------------------------------------
%           Initializing
for tb=1:1000
    
    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        mincf=min(F);
        maxcf=max(F);
        normf=(F-mincf)/(maxcf-mincf);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        ansc=evalfis([i*ones(nSW,1) normf],nbn);
        acc_coef1(acc_coef1~=0)=ansc(:,1);
        acc_coef2(acc_coef2~=0)=ansc(:,2);
        Vel=inertia_weight*Vel+acc_coef1*(rand(nSW,nPam).*(Pbest-X))+acc_coef2*(rand(nSW,nPam).*(ones(nSW,1)*Gbest-X));
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    fuzpsoGriewangk8r3(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
end
fprintf(' Ýter=%d ObjVal=%g\n',var(fuzpsoGriewangk8r3),mean(fuzpsoGriewangk8r3))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('fuzpsoGriewangk8r3','fuzpsoGriewangk8r3')





clear all;
close all;
%------------- PSO code
%% Parameter
nSW=50;  % number of particles
nPam=8; % dimension; number of design variable
nIter_max=100;  % maximum of iterations
min_fit_criteria=0; % acceptable minimum value of fitness
inertia_weight=0.72;
acc_coef1=eye(nSW);
acc_coef2=eye(nSW);
ub=600;
lb=-600;
ubt=600;
lbt=-600;

%---------------------------------------------------
%           fuzzy Parameter
nbn=newfis('fpso');
nbn=addvar(nbn,'input','iter',[0 nIter_max]);
nbn=addmf(nbn,'input',1,'start','gaussmf',[.25*nIter_max 0]);
nbn=addmf(nbn,'input',1,'mid2','gaussmf',[.25*nIter_max .5*nIter_max]);
nbn=addmf(nbn,'input',1,'end','gaussmf',[.25*nIter_max nIter_max]);

nbn=addvar(nbn,'input','cfval',[0 1]);
nbn=addmf(nbn,'input',2,'bad','gaussmf',[.25 0]);
nbn=addmf(nbn,'input',2,'normal','gaussmf',[.25 .5]);
nbn=addmf(nbn,'input',2,'excellent','gaussmf',[.25 1]);

nbn=addvar(nbn,'output','c1',[0 3]);
nbn=addmf(nbn,'output',1,'lowc1','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',1,'medium1c1','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',1,'medium2c1','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',1,'medium3c1','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',1,'highc1','trimf',[2 2.5 3]);

nbn=addvar(nbn,'output','c2',[0 3]);
nbn=addmf(nbn,'output',2,'lowc2','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',2,'medium1c2','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',2,'medium2c2','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',2,'medium3c2','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',2,'highc2','trimf',[2 2.5 3]);

% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 2 4 1 2
%     1 3 1 5 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 5 1 1 2
%     3 2 4 2 1 2
%     3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 4 2 1 2
%     1 3 5 1 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 1 5 1 2
%     3 2 2 4 1 2
%     3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 1 5 1 2
%     1 2 2 4 1 2
%     1 3 3 3 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 3 3 1 2
%     3 2 4 2 1 2
%     3 3 5 1 1 2];
ruleList=[ ...
    1 1 5 1 1 2
    1 2 4 2 1 2
    1 3 3 3 1 2
    2 1 4 2 1 2
    2 2 3 3 1 2
    2 3 2 4 1 2
    3 1 3 3 1 2
    3 2 2 4 1 2
    3 3 1 5 1 2];

nbn=addrule(nbn,ruleList);

%  -----------------------------------------------------
%           Initializing
for tb=1:1000
    
    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        mincf=min(F);
        maxcf=max(F);
        normf=(F-mincf)/(maxcf-mincf);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        ansc=evalfis([i*ones(nSW,1) normf],nbn);
        acc_coef1(acc_coef1~=0)=ansc(:,1);
        acc_coef2(acc_coef2~=0)=ansc(:,2);
        Vel=inertia_weight*Vel+acc_coef1*(rand(nSW,nPam).*(Pbest-X))+acc_coef2*(rand(nSW,nPam).*(ones(nSW,1)*Gbest-X));
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    fuzpsoGriewangk8r4(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
end
fprintf(' Ýter=%d ObjVal=%g\n',var(fuzpsoGriewangk8r4),mean(fuzpsoGriewangk8r4))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('fuzpsoGriewangk8r4','fuzpsoGriewangk8r4')





clear all;
close all;
%------------- PSO code
    %% Parameter
    nSW=50;  % number of particles
    nPam=20; % dimension; number of design variable
    nIter_max=100;  % maximum of iterations
    min_fit_criteria=0; % acceptable minimum value of fitness
    inertia_weight=0.72;
    acc_coef1=2;
    acc_coef2=2;
    ub=600;
    lb=-600;
    ubt=600;
    lbt=-600;
    %  -----------------------------------------------------
    %           Initializing
        for tb=1:1000

    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    %---------------------------------------------------
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        Vel=inertia_weight*Vel+acc_coef1*rand(nSW,nPam).*(Pbest-X)+acc_coef2*rand(nSW,nPam).*(ones(nSW,1)*Gbest-X);
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    psoGriewangk20(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
        end
fprintf(' Ýter=%d ObjVal=%g\n',var(psoGriewangk20),mean(psoGriewangk20))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('psoGriewangk20','psoGriewangk20')





clear all;
close all;
%------------- PSO code
%% Parameter
nSW=50;  % number of particles
nPam=20; % dimension; number of design variable
nIter_max=100;  % maximum of iterations
min_fit_criteria=0; % acceptable minimum value of fitness
inertia_weight=0.72;
acc_coef1=eye(nSW);
acc_coef2=eye(nSW);
ub=600;
lb=-600;
ubt=600;
lbt=-600;

%---------------------------------------------------
%           fuzzy Parameter
nbn=newfis('fpso');
nbn=addvar(nbn,'input','iter',[0 nIter_max]);
nbn=addmf(nbn,'input',1,'start','gaussmf',[.25*nIter_max 0]);
nbn=addmf(nbn,'input',1,'mid2','gaussmf',[.25*nIter_max .5*nIter_max]);
nbn=addmf(nbn,'input',1,'end','gaussmf',[.25*nIter_max nIter_max]);

nbn=addvar(nbn,'input','cfval',[0 1]);
nbn=addmf(nbn,'input',2,'bad','gaussmf',[.25 0]);
nbn=addmf(nbn,'input',2,'normal','gaussmf',[.25 .5]);
nbn=addmf(nbn,'input',2,'excellent','gaussmf',[.25 1]);

nbn=addvar(nbn,'output','c1',[0 3]);
nbn=addmf(nbn,'output',1,'lowc1','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',1,'medium1c1','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',1,'medium2c1','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',1,'medium3c1','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',1,'highc1','trimf',[2 2.5 3]);

nbn=addvar(nbn,'output','c2',[0 3]);
nbn=addmf(nbn,'output',2,'lowc2','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',2,'medium1c2','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',2,'medium2c2','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',2,'medium3c2','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',2,'highc2','trimf',[2 2.5 3]);

ruleList=[ ...
    1 1 3 3 1 2
    1 2 2 4 1 2
    1 3 1 5 1 2
    2 1 4 2 1 2
    2 2 3 3 1 2
    2 3 2 4 1 2
    3 1 5 1 1 2
    3 2 4 2 1 2
    3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 4 2 1 2
%     1 3 5 1 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 1 5 1 2
%     3 2 2 4 1 2
%     3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 1 5 1 2
%     1 2 2 4 1 2
%     1 3 3 3 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 3 3 1 2
%     3 2 4 2 1 2
%     3 3 5 1 1 2];
% ruleList=[ ...
%     1 1 5 1 1 2
%     1 2 4 2 1 2
%     1 3 3 3 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 3 3 1 2
%     3 2 2 4 1 2
%     3 3 1 5 1 2];

nbn=addrule(nbn,ruleList);

%  -----------------------------------------------------
%           Initializing
for tb=1:1000
    
    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        mincf=min(F);
        maxcf=max(F);
        normf=(F-mincf)/(maxcf-mincf);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        ansc=evalfis([i*ones(nSW,1) normf],nbn);
        acc_coef1(acc_coef1~=0)=ansc(:,1);
        acc_coef2(acc_coef2~=0)=ansc(:,2);
        Vel=inertia_weight*Vel+acc_coef1*(rand(nSW,nPam).*(Pbest-X))+acc_coef2*(rand(nSW,nPam).*(ones(nSW,1)*Gbest-X));
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    fuzpsoGriewangk20r1(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
end
fprintf(' Ýter=%d ObjVal=%g\n',var(fuzpsoGriewangk20r1),mean(fuzpsoGriewangk20r1))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('fuzpsoGriewangk20r1','fuzpsoGriewangk20r1')






clear all;
close all;
%------------- PSO code
%% Parameter
nSW=50;  % number of particles
nPam=20; % dimension; number of design variable
nIter_max=100;  % maximum of iterations
min_fit_criteria=0; % acceptable minimum value of fitness
inertia_weight=0.72;
acc_coef1=eye(nSW);
acc_coef2=eye(nSW);
ub=600;
lb=-600;
ubt=600;
lbt=-600;

%---------------------------------------------------
%           fuzzy Parameter
nbn=newfis('fpso');
nbn=addvar(nbn,'input','iter',[0 nIter_max]);
nbn=addmf(nbn,'input',1,'start','gaussmf',[.25*nIter_max 0]);
nbn=addmf(nbn,'input',1,'mid2','gaussmf',[.25*nIter_max .5*nIter_max]);
nbn=addmf(nbn,'input',1,'end','gaussmf',[.25*nIter_max nIter_max]);

nbn=addvar(nbn,'input','cfval',[0 1]);
nbn=addmf(nbn,'input',2,'bad','gaussmf',[.25 0]);
nbn=addmf(nbn,'input',2,'normal','gaussmf',[.25 .5]);
nbn=addmf(nbn,'input',2,'excellent','gaussmf',[.25 1]);

nbn=addvar(nbn,'output','c1',[0 3]);
nbn=addmf(nbn,'output',1,'lowc1','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',1,'medium1c1','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',1,'medium2c1','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',1,'medium3c1','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',1,'highc1','trimf',[2 2.5 3]);

nbn=addvar(nbn,'output','c2',[0 3]);
nbn=addmf(nbn,'output',2,'lowc2','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',2,'medium1c2','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',2,'medium2c2','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',2,'medium3c2','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',2,'highc2','trimf',[2 2.5 3]);

% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 2 4 1 2
%     1 3 1 5 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 5 1 1 2
%     3 2 4 2 1 2
%     3 3 3 3 1 2];
ruleList=[ ...
    1 1 3 3 1 2
    1 2 4 2 1 2
    1 3 5 1 1 2
    2 1 2 4 1 2
    2 2 3 3 1 2
    2 3 4 2 1 2
    3 1 1 5 1 2
    3 2 2 4 1 2
    3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 1 5 1 2
%     1 2 2 4 1 2
%     1 3 3 3 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 3 3 1 2
%     3 2 4 2 1 2
%     3 3 5 1 1 2];
% ruleList=[ ...
%     1 1 5 1 1 2
%     1 2 4 2 1 2
%     1 3 3 3 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 3 3 1 2
%     3 2 2 4 1 2
%     3 3 1 5 1 2];

nbn=addrule(nbn,ruleList);

%  -----------------------------------------------------
%           Initializing
for tb=1:1000
    
    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        mincf=min(F);
        maxcf=max(F);
        normf=(F-mincf)/(maxcf-mincf);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        ansc=evalfis([i*ones(nSW,1) normf],nbn);
        acc_coef1(acc_coef1~=0)=ansc(:,1);
        acc_coef2(acc_coef2~=0)=ansc(:,2);
        Vel=inertia_weight*Vel+acc_coef1*(rand(nSW,nPam).*(Pbest-X))+acc_coef2*(rand(nSW,nPam).*(ones(nSW,1)*Gbest-X));
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    fuzpsoGriewangk20r2(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
end
fprintf(' Ýter=%d ObjVal=%g\n',var(fuzpsoGriewangk20r2),mean(fuzpsoGriewangk20r2))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('fuzpsoGriewangk20r2','fuzpsoGriewangk20r2')






clear all;
close all;
%------------- PSO code
%% Parameter
nSW=50;  % number of particles
nPam=20; % dimension; number of design variable
nIter_max=100;  % maximum of iterations
min_fit_criteria=0; % acceptable minimum value of fitness
inertia_weight=0.72;
acc_coef1=eye(nSW);
acc_coef2=eye(nSW);
ub=600;
lb=-600;
ubt=600;
lbt=-600;

%---------------------------------------------------
%           fuzzy Parameter
nbn=newfis('fpso');
nbn=addvar(nbn,'input','iter',[0 nIter_max]);
nbn=addmf(nbn,'input',1,'start','gaussmf',[.25*nIter_max 0]);
nbn=addmf(nbn,'input',1,'mid2','gaussmf',[.25*nIter_max .5*nIter_max]);
nbn=addmf(nbn,'input',1,'end','gaussmf',[.25*nIter_max nIter_max]);

nbn=addvar(nbn,'input','cfval',[0 1]);
nbn=addmf(nbn,'input',2,'bad','gaussmf',[.25 0]);
nbn=addmf(nbn,'input',2,'normal','gaussmf',[.25 .5]);
nbn=addmf(nbn,'input',2,'excellent','gaussmf',[.25 1]);

nbn=addvar(nbn,'output','c1',[0 3]);
nbn=addmf(nbn,'output',1,'lowc1','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',1,'medium1c1','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',1,'medium2c1','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',1,'medium3c1','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',1,'highc1','trimf',[2 2.5 3]);

nbn=addvar(nbn,'output','c2',[0 3]);
nbn=addmf(nbn,'output',2,'lowc2','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',2,'medium1c2','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',2,'medium2c2','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',2,'medium3c2','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',2,'highc2','trimf',[2 2.5 3]);

% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 2 4 1 2
%     1 3 1 5 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 5 1 1 2
%     3 2 4 2 1 2
%     3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 4 2 1 2
%     1 3 5 1 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 1 5 1 2
%     3 2 2 4 1 2
%     3 3 3 3 1 2];
ruleList=[ ...
    1 1 1 5 1 2
    1 2 2 4 1 2
    1 3 3 3 1 2
    2 1 2 4 1 2
    2 2 3 3 1 2
    2 3 4 2 1 2
    3 1 3 3 1 2
    3 2 4 2 1 2
    3 3 5 1 1 2];
% ruleList=[ ...
%     1 1 5 1 1 2
%     1 2 4 2 1 2
%     1 3 3 3 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 3 3 1 2
%     3 2 2 4 1 2
%     3 3 1 5 1 2];

nbn=addrule(nbn,ruleList);

%  -----------------------------------------------------
%           Initializing
for tb=1:1000
    
    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        mincf=min(F);
        maxcf=max(F);
        normf=(F-mincf)/(maxcf-mincf);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        ansc=evalfis([i*ones(nSW,1) normf],nbn);
        acc_coef1(acc_coef1~=0)=ansc(:,1);
        acc_coef2(acc_coef2~=0)=ansc(:,2);
        Vel=inertia_weight*Vel+acc_coef1*(rand(nSW,nPam).*(Pbest-X))+acc_coef2*(rand(nSW,nPam).*(ones(nSW,1)*Gbest-X));
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    fuzpsoGriewangk20r3(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
end
fprintf(' Ýter=%d ObjVal=%g\n',var(fuzpsoGriewangk20r3),mean(fuzpsoGriewangk20r3))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('fuzpsoGriewangk20r3','fuzpsoGriewangk20r3')





clear all;
close all;
%------------- PSO code
%% Parameter
nSW=50;  % number of particles
nPam=20; % dimension; number of design variable
nIter_max=100;  % maximum of iterations
min_fit_criteria=0; % acceptable minimum value of fitness
inertia_weight=0.72;
acc_coef1=eye(nSW);
acc_coef2=eye(nSW);
ub=600;
lb=-600;
ubt=600;
lbt=-600;

%---------------------------------------------------
%           fuzzy Parameter
nbn=newfis('fpso');
nbn=addvar(nbn,'input','iter',[0 nIter_max]);
nbn=addmf(nbn,'input',1,'start','gaussmf',[.25*nIter_max 0]);
nbn=addmf(nbn,'input',1,'mid2','gaussmf',[.25*nIter_max .5*nIter_max]);
nbn=addmf(nbn,'input',1,'end','gaussmf',[.25*nIter_max nIter_max]);

nbn=addvar(nbn,'input','cfval',[0 1]);
nbn=addmf(nbn,'input',2,'bad','gaussmf',[.25 0]);
nbn=addmf(nbn,'input',2,'normal','gaussmf',[.25 .5]);
nbn=addmf(nbn,'input',2,'excellent','gaussmf',[.25 1]);

nbn=addvar(nbn,'output','c1',[0 3]);
nbn=addmf(nbn,'output',1,'lowc1','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',1,'medium1c1','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',1,'medium2c1','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',1,'medium3c1','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',1,'highc1','trimf',[2 2.5 3]);

nbn=addvar(nbn,'output','c2',[0 3]);
nbn=addmf(nbn,'output',2,'lowc2','trimf',[0 .5 1]);
nbn=addmf(nbn,'output',2,'medium1c2','trimf',[.5 1 1.5]);
nbn=addmf(nbn,'output',2,'medium2c2','trimf',[1 1.5 2]);
nbn=addmf(nbn,'output',2,'medium3c2','trimf',[1.5 2 2.5]);
nbn=addmf(nbn,'output',2,'highc2','trimf',[2 2.5 3]);

% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 2 4 1 2
%     1 3 1 5 1 2
%     2 1 4 2 1 2
%     2 2 3 3 1 2
%     2 3 2 4 1 2
%     3 1 5 1 1 2
%     3 2 4 2 1 2
%     3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 3 3 1 2
%     1 2 4 2 1 2
%     1 3 5 1 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 1 5 1 2
%     3 2 2 4 1 2
%     3 3 3 3 1 2];
% ruleList=[ ...
%     1 1 1 5 1 2
%     1 2 2 4 1 2
%     1 3 3 3 1 2
%     2 1 2 4 1 2
%     2 2 3 3 1 2
%     2 3 4 2 1 2
%     3 1 3 3 1 2
%     3 2 4 2 1 2
%     3 3 5 1 1 2];
ruleList=[ ...
    1 1 5 1 1 2
    1 2 4 2 1 2
    1 3 3 3 1 2
    2 1 4 2 1 2
    2 2 3 3 1 2
    2 3 2 4 1 2
    3 1 3 3 1 2
    3 2 2 4 1 2
    3 3 1 5 1 2];

nbn=addrule(nbn,ruleList);

%  -----------------------------------------------------
%           Initializing
for tb=1:1000
    
    Vel=zeros(nSW,nPam);
    X=rand(nSW,nPam)*(ubt-lbt)+lbt;
    Pbest=X;
    Gbest=X(1,:);
    F0=1e10*ones(nSW,1);
    min_F=1e10;
    i=0;
    %% Evolutionary Loop
    while (i<nIter_max)&&(min_F>min_fit_criteria)
        i=i+1;
        F=fit_fun_Griewangk(X);
        mincf=min(F);
        maxcf=max(F);
        normf=(F-mincf)/(maxcf-mincf);
        index=find((F-F0)<0);
        Pbest(index,:)=X(index,:);
        F0(index)=F(index);
        [min_F,ind_Gbest]=min(F0);
        Gbest=Pbest(ind_Gbest,:);
        ansc=evalfis([i*ones(nSW,1) normf],nbn);
        acc_coef1(acc_coef1~=0)=ansc(:,1);
        acc_coef2(acc_coef2~=0)=ansc(:,2);
        Vel=inertia_weight*Vel+acc_coef1*(rand(nSW,nPam).*(Pbest-X))+acc_coef2*(rand(nSW,nPam).*(ones(nSW,1)*Gbest-X));
        Vel=min(Vel,30);
        Vel=max(Vel,-30);
        X=X+Vel;
        X=max(X,-600); X=min(X,600);
    end
    fuzpsoGriewangk20r4(tb)=min_F;
%     fprintf('Ýter=%d ObjVal=%g\n',tb,min_F);
end
fprintf(' Ýter=%d ObjVal=%g\n',var(fuzpsoGriewangk20r4),mean(fuzpsoGriewangk20r4))
% mean(fuzpsoAckley)
% var(fuzpsoAckley)
save('fuzpsoGriewangk20r4','fuzpsoGriewangk20r4')
