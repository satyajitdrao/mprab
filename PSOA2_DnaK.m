% function [Solution,x] = PSOA2_DnaK
% parameter selection for model with explicit DnaK. Calls @errA2_DnaK
% and @modelA2_DnaK
% clear all;
% global lb ub
         % 1     2     3     4     5      6      7
         % kpdeg kmdeg ktlnA ktlnR kbtpn1 kbtpn2 kbtpn3 
lb(1:7) = [-4.8	 -2.92 -2    1     -6	  -6     -6.2]; 
ub(1:7) = [-4.4	 -2.92 -1.4  2     -5.8	  -5.7   -6];
        %  8    9   10   11   12   13    14   15   16   17   18    19 20   21   22
        %  f1   f2  f3   K1   K2   ktlnD SDS  kdsf K5b  kdsb krdeg kb kd   kbdf kbdb
lb(8:22) = [0.4 0.5 1   -3.5 -1.8  0     1    -1   -2   -3.8 -4.6  -1 -3.8 -1   -3.8];
ub(8:22) = [0.7 1   1.3 -2   -0.5  2     1.7   0   1    -3.0 -3.7   1 -1.5  0   -3];
            %23   24  25  26   27   28   29   30    31       32   33     34  35
            %k2   km2 k3  km3  k4   k5   km5  k6   factor2  kfit1 kfit2 n1  n2
lb(23:35) =[-2.5  -5 -1.4 -3.7 -2  -1.5 -3.7 -2.5  1.2      -2    -2    0   0];
ub(23:35) =[ 0    -3  0   0.5  0.5  0.3 -2   -1	   2        -1    -1    1   1];

% % % % % % Optimizing MprA --> mprA/sigE
            %  f1b f3b
lb([36 37]) = [0   0];
ub([36 37]) = [0.7 1.3];

% % % % % % changed: no krdeg; set x(18) = -Inf in errA2_DnaK;
% lb(31) = -1; ub(31) = -1;
% lb(18) = -1; ub(18) = -1;

n_param=length(lb);
nsol=50;
Solution=ones(nsol+2,n_param+1)*10000;
Solution(nsol+1,:)=cat(2,lb, 10);
Solution(nsol+2,:)=cat(2,ub, 10);
for j=1:nsol
    options = optimoptions('particleswarm','SwarmSize',50,'Display','iter','PlotFcn',@pswplotbestf);
    options.MaxIterations = 80; options.FunctionTolerance = 5e-3;
   [X, fval] = particleswarm(@errA2_DnaK,n_param,lb,ub,options);
   Solution(j,:) = [X fval];
   save('09202020','Solution')
end



% lb(2) = -2.94; ub(2) = -2.90; ub(16) = 1.2;
% options2 = optimoptions('fmincon','PlotFcn',@optimplotfval,'Display','iter','OptimalityTolerance',1e-3);
% [s, fval] = fmincon(@errA2_DnaK,x,[],[],[],[],lb,ub,[],options2)