function [Solution,x] = PSOA2

clear all;
global lb ub
% parameters = ['kpdeg kmdeg ktlnA ktlnB ktlnE  ktlnR kbtpn1 kbtpn2 kbtpn3' ...
%     ' f1 f2 f3 K1 K2 K3 f5a f5b K5a K5b f6 K4 kb4 krdeg kb kd k1 km1' ...
%     ' k2 km2 k3 km3 k4 k5 km5 k6 kb2 kd2 krp krd factor1 factor2'];
     % 1     2     3     4     5      6      7
     % kpdeg kmdeg ktlnA ktlnR kbtpn1 kbtpn2 kbtpn3 
lb1 = [-4.8	 -2.92 -2.0  0     -6	  -6     -6.3]; 
ub1 = [-4.1	 -2.92 -1.4  1.7   -5.7	  -5.7   -5.85];
    % 8    9   10  11   12   13   14    15 16   
    % f1   f2  f3  K1   K2   K5b  krdeg kb kd   
lb2 = [0.4 0.5 0.7 -3.5 -1.8 0.5  -5    -2 -3.8 ];
ub2 = [0.7 1   1.4 -2   -0.5 1.4  -3     1 -1.5 ];
     %17   18  19   20  21 22   23  24   25   26    27      28      29    30    31  32
     %k10  km1 k2   km2 k3 km3  k4  k5   km5  k6    factor1 factor2 kfit1 kfit2 n1  n2
lb3 =[-3.7 -3  -2.5 -5  -2 -3   -2  -1.5 -3.7 -2.5  1.7	    0.3     -3    -3    0   0];
ub3 =[-1.3 -1  -1   -2  0  0.5  0.5 0.3  -2   -1	3	    3        0     0    1.3 1.3];


lb = [lb1 lb2 lb3];
ub = [ub1 ub2 ub3];

% % % remove K5b
lb(13) = 0; ub(13) = 0;
% % % 

n_param=length(lb);
nrep=50;
nsol=min(nrep,20);
Solution=ones(nsol+2,n_param+1)*10000;
Solution(nsol+1,:)=cat(2,lb, 10);
Solution(nsol+2,:)=cat(2,ub, 10);
for j=1:nrep

    options = optimoptions('particleswarm','SwarmSize',50,'Display','iter'); %'PlotFcn',@pswplotbestf,
    options.MaxIterations = 100; options.FunctionTolerance = 1e-3;
   [X, fval] = particleswarm(@errA2,n_param,lb,ub,options);
   Solution(j,:) = [X fval];
   Solution=sortrows(Solution,n_param+1);
end
   save('0531','Solution')
end
% x=Solution(1,:);
% % x = [-4.11346013669348,-2.92000000000000,-1.95561030065855,1.24066089107222,-5.84376121525913,-5.76493491512478,-6.08961890770863,0.459872333516754,0.819266311124835,1.02215064476271,-2.80542296088934,-1.45245454893285,0,-4.19735454831023,0.987010356726607,-3.06929457796160,-3.36277653879768,-1.40932539802979,-2.11085831253633,-4.44150546900535,-1.99960286213257,0.354029025709110,-1.35127860742104,-1.15593685402913,-2.74980176753139,-2.35451410935175,2.09930652342334,2.50815069616319,-2.00117297522030,-1.31190164081655,0.487959937441361,0.813947359797499,0.266563189052063];
% % options2 = optimoptions('fmincon','PlotFcn',@optimplotfval,'Display','iter');
% % fmincon(@errA2,x,[],[],[],[],lb,ub,[],options2)