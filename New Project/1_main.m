%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Dan Bu, Date: April 21th, 2015
% Obj: solve the semidefinite relaxation of the fractional programming
% model in production-inventory integrated system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%% vary K
objK = [];
num = 0;
for i=100:20:140%2080
    num = num + 1;
    K = [300,i];
    objK(num,:) = varyK(K);
end
% obj_fr_value/obj_fs_value
ratioK = resultK(:,2)./resultK(:,1);
%% vary h
objH = [];
num = 0;
for i=11:110
    num = num + 1;
    h = [10,i];
    [a,b] = varyH(h);
    objH(num,:) = [a,b];
end
ratioH = resultH(:,2)./resultH(:,1);
%% display
plot(ratioK)
plot(ratioH)

%% outputs
csvwrite('ratioK.csv',ratioK)
csvwrite('ratioH.csv',ratioH)