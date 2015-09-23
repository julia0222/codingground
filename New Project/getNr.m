%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Dan Bu, Date: April 21th, 2015
% Obj: solve the semidefinite relaxation of the fractional programming
% model in production-inventory integrated system
% fixed ratio shipment policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[n] = getNr(K,h,D,mu,L,Q)
% K = [1000, 200];  % add '.' to indicate that it's a real number
% h = [5, 10];
% D = [10,10];
% mu = [30., 40.];
% p = [0.3, 0.7];   % probability distribution of mu
% L = [50., 20.];    % perishability length
% Q = [50,500];        % number of batches per production cycle
% Q = tmpQ

% calculate n
k = mu./D;
K1 = K(2);
C1 = h(1)*Q.^2.*(k-1)./mu;
C2 = (h(2)-h(1)).*Q.^2.*(k-1)./(2*D.*(k+1));
M = K1./(log(k).*(C1+2*C2));
x = zeros(length(mu),2);
a = M;
b = -(2*M+1);
c = M;
d = sqrt(b.^2 - 4*a.*c);
x(:,1) = ( -b + d ) ./ (2*a);
x(:,2) = ( -b - d ) ./ (2*a);

if length(x(x>1)) ~= length(mu)
    error('double check');
end

n = log(x(x>1))'./log(k);
nfloor = floor(n);
nceil = ceil(n);
v1 = f(nfloor,k,K1,C1,C2);
v2 = f(nceil,k,K1,C1,C2);

% calculate lower bound for n
tmp1 = L(1)*k.*mu./((k-1).*Q);
tmp2 = L(2)*k.*D./((k-1).*Q);
UB = min(tmp1,tmp2);
LB = log(UB./(UB-1))./log(k);

for i=1:length(mu)
    if nfloor(i) < LB(i)
       n(i) = ceil(LB(i));
    elseif v1(i) < v2(i)
       n(i) = nfloor(i);
    else  
       n(i) = nceil(i);
    end
end


