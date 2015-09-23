%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Dan Bu, Date: April 21th, 2015
% Obj: solve the semidefinite relaxation of the fractional programming
% model in production-inventory integrated system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%% Parameters
num = 100;
K = [300,100];  % 
h = [10,20];
D = [5,5,5];
mu = [40, 30, 20];
p = [0.3, 0.2, 0.5];   % probability distribution of mu
L = [20,50];    % perishability length
epsi = 0.01;

%% fixed size shipment optimization - no perishability
idx_fs = 1;
fsN = ones(1,length(mu));
fsQ = zeros(1,length(mu));
while idx_fs <= num
    W = getQ(K,h,D,mu,p,L,fsN(end,:)); 
    %[W0, valueH0] = getQ0(K,h,D,mu,p,L,N(end,:));
    if sum(eigs(W) < 0.0001) == length(mu)
        [Vec, Val] = eigs(W, 1);
        tmpQ = (Vec(1:length(mu))/Vec(end))';
    else 
        error('not feasible');  
        break;
    end
    tmpN = getN(K,h,D,mu,L,tmpQ);
    
    % termination condition
    if sum(abs(tmpQ-fsQ(end,:)) < epsi) == length(mu)
       break;
    else
       fsQ = [fsQ;tmpQ];
       fsN = [fsN;tmpN];
    end
    idx_fs = idx_fs + 1;
end
obj_fs_value = obj_fs(K,h,mu,D,p,fsQ(end,:),fsN(end,:))

%% fixed ratio shipment optimization - no perishability
idx_fr = 1;
k = mu./D;
frN = ones(1,length(mu));
frQ = zeros(1,length(mu));
while idx_fr <= num
    [W, U] = getQr(K,h,D,mu,p,L,frN(end,:));
    if sum(eigs(W) < 0.0001) == length(mu)
        [Vec, Val] = eigs(W, 1);
        tmpQ = (Vec(1:length(mu))/Vec(end))';
    else 
        error('fixed ratio not feasible');
        break;
    end
    tmpN = getNr(K,h,D,mu,L,tmpQ); 
    
    % termination condition
    if sum(abs(tmpQ-frQ(end,:)) < epsi) == length(mu)
       break;
    else
       frQ = [frQ;tmpQ];
       frN = [frN;tmpN];
    end
    idx_fr = idx_fr + 1;
end
obj_fr_value = obj_fr(K,h,mu,D,p,frQ(end,:),frN(end,:))

%% display




