%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Dan Bu, Date: April 21th, 2015
% Obj: solve the semidefinite relaxation of the fractional programming
% model in production-inventory integrated system
% fixed ratio shipment policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[W, U] = getQr(K,h,D,mu,p,L,n)
k = mu./D;
% calculate U -- upper bound of Q(i)
U = zeros(1,length(n));
tmp1 = L(1)*mu.*(k.^n-1)./(k.^(n-1).*(k-1));
tmp2 = L(2)*D.*(k.^n-1)./(k.^(n-1).*(k-1));
tmp3 = L(2)*k.*D./(k-1);   % got from "optimize n given Q" to gurantee the feasibility of n
U = min(min(tmp1,tmp2),tmp3);
clear tmp1 tmp2;

% calculate A
A_tild = diag(p.*(h(1)/2.*(1./D-1./mu) + h(1)*(k-1)./(mu.*(k.^n-1)) + (h(2)-h(1))./(2*D).*((k-1).*(k.^n+1))./((k+1).*(k.^n-1))));
b = zeros(length(n),1);
c = K(1) + K(2)*sum(n.*p);
A = [A_tild, b; b',c];

% calculate B
B_tild = zeros(length(n));
b = 1/2*(p./D)';
c = 0;
B = [B_tild, b; b',c];

% calculate C
C = zeros(length(n)+1 ,length(n)+1, length(n));
for i=1:length(n)
    vec_a = zeros(length(n),1)';
    vec_a(i) = 1;
    C_tild = vec_a'*vec_a;
    b = -1/2*U(i)*vec_a;
    c = 0;
    C(:,:,i) = [C_tild, b'; b,c];   
end

%% Optimization
cvx_begin sdp quiet
% define variables
variable W(length(n)+1 ,length(n)+1) symmetric

% optimization
minimize(trace(A*W));

% add constraints
subject to
    trace(B*W) >= 1;
    trace(B*W) <= 1;
    for i=1:length(n)
        trace(C(:,:,i)*W) <= 0;
    end
    W == semidefinite(length(n)+1);
cvx_end

value = trace(A*W);
end