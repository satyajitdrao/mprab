function out=errA2(x)
% global lb ub
global ktlnB k1 km1 k2 km2 k3 km3 k4 k5 km5 k6
global kpdeg kmdeg ktlnA ktlnE ktlnR kbtpn1 kbtpn2 kbtpn3 f1 f2 f3 K1 K2 K3 kb kd krdeg
global f5b K5b

% if min(logical(lb<=x & ub>=x)==1)
%      A  Ap   E   R  ER  mAB   mE    Bs  B   Bp  kpa   pap   sigB
X0 = zeros(1,13); %[0.1 0.1 0.1 0.1 0.1 0.001 0.001 0.1 0.1 0.1 0.001 0.001 0.001];

% translation, transcription
kpdeg= 10^x(1);
kmdeg = 10^x(2);
ktlnA = 10^x(3); ktlnA0 = ktlnA;
ktlnB = ktlnA0/4.5;
ktlnE = ktlnA0/5; ktlnE0 = ktlnE;
kbtpn1 = 10^x(5);
kbtpn2 = 10^x(6);
kbtpn3 = 10^x(7);
ktlnR = kbtpn3*ktlnE*10^x(4)/kmdeg;
f1 = 10^x(8);
f2  = 10^x(9);
f3  = 10^x(10);
K1 = 10^x(11);
K2 = 10^x(12);
K3 = K1;
    % MprA-DNA binding and transcription initiation
f5b = f3;
K5b = Inf; %10^x(13);
% post-translational
krdeg0 = 10^x(14);
kb = 10^x(15);
kd = 10^x(16);
k10 = 10^x(17);
km1 = 10^x(18);
k2 = 10^x(19);
km2 = 10^x(20);
k3= 10^x(21);
km3= 10^x(22);
k4= 10^x(23);
k5 = 10^x(24);
km5 = 10^x(25);
k6 = 10^x(26);
factor1 = 10^x(27);
factor2 = 10^x(28);
Kfit1 = 10^x(29); Kfit2 = 10^x(30);
n1 = 10^x(31); n2= 10^x(32);

krdegf = @(SDS) krdeg0 + krdeg0*factor2*(SDS.^n2./(SDS.^n2 + Kfit2^n2));
k1f = @(SDS) k10+k10*factor1*(SDS.^n1./(SDS.^n1 + Kfit1^n1));
%% Dynamics 0,0.03
    tdif = 60;
tspan = 0:tdif:2*3600;
% pre-stress steady state
k1=k10;
krdeg=krdeg0;
[~, X]=ode15s(@modelA2, [0 10*3600], X0); % 0 ss
% bacteriostatic stress - 2 hrs @ 0.03%
k1 = k1f(0.03);
krdeg = krdegf(0.03);
[~,Y] = ode15s(@modelA2, tspan, X(end,:)); % 0 --> 0.03
k1=k10;
krdeg=krdeg0;
[~, X2]=ode15s(@modelA2, tspan, Y(end,:)); % 0.03 --> 0

SDS = [0,0.01,0.015,0.02,0.025,0.03];
UP = zeros(length(SDS),length(X0));Down = zeros(length(SDS),length(X0));
for i = 1:length(SDS)
k1 = k1f(SDS(i));
krdeg = krdegf(SDS(i));
[~,U] = ode15s(@modelA2, tspan, X(end,:));
[~,D] = ode15s(@modelA2, tspan, Y(end,:));
UP(i,:) = U(end,:);
Down(i,:) = D(end,:);
end
UP(1,:) = X(end,:);
Down(1,:) = X2(end,:);
%% data
err = zeros(1,12);
idx_u = 1+[0 15 30 60 120]*60/tdif;
idx_d = 1+[0 30 60 120]*60/tdif;

sigE0 = [0.007816802	0.008776893	0.010462324	0.008700823	0.010163107	0.006815801];
sigB0 = [0.009949661	0.009981423	0.014413074	0.005179989	0.005366144	0.003649808];
mprA0 = [0.000892265	0.001087373	0.001027269	0.000502307	0.000491452	0.000332753];
% dynamic data
sigE_0_30 = [0	mean(sigE0) mean(sigE0) mean(sigE0)	
15	0.047338157	0.036133309	0.040389475			
30	0.099011552	0.107391634	0.095859325			
60	0.091994471	0.110281237	0.0907232			
120	0.087854034	0.096661217	0.092778504];			
						

sigE_30_0 = [0	0.087854034	0.096661217	0.092778504			
30	0.026582011	0.034363073	0.023439055			
60	0.014266473	0.024790752	0.025814532			
120	0.021392016	0.021084555	0.018204148];
			
						
% sigB						
sigB_0_30 =	[0 mean(sigB0) mean(sigB0) mean(sigB0)	
15	0.025735146	0.020825684	0.023756414			
30	0.126505011	0.119252291	0.103545824			
60	0.260074617	0.229640044	0.23244971			
120	0.153855923	0.195462021	0.190877581];			
						
						
sigB_30_0 = [0	0.153855923	0.195462021	0.190877581			
30	0.029783596	0.044297444	0.028836135			
60	0.006975988	0.005821612	0.007264636			
120	0.007641807	0.008628851	0.007792603];
						
% mprA						
mprA_0_30 = [0	mean(mprA0) mean(mprA0) mean(mprA0)
15	0.003558258	0.002549538	0.003054443			
30	0.004887063	0.004685978	0.004603986			
60	0.006798646	0.006143091	0.004752837			
120	0.006175092	0.006496226	0.006628581];			
						
						
mprA_30_0 = [0	0.006175092	0.006496226	0.006628581			
30	0.002863715	0.003970514	0.002902044			
60	0.000530308	0.000707476	0.000915335			
120	0.000563387	0.000708367	0.000691986];		

% dose response data
sigE_sds_inc3 = [0	mean(sigE0) mean(sigE0) mean(sigE0)	
0.01	0.008631604	0.009595789	0.008289541			
0.015	0.01222866	0.00858139	0.012360231			
0.02	0.015974817	0.029355593	0.026730883			
0.025	0.062927849	0.060524653	0.058466244			
0.03	0.087854034	0.096661217	0.092778504];			
sigE_sds_dec3 = [0.03	0.0898023	0.090497034	0.079647766			
0.025	0.094679283	0.071230185	0.088266159			
0.02	0.095120458	0.069976793	0.081862344			
0.015	0.023464665	0.022023889	0.025056699			
0.01	0.023350733	0.021475609	0.019307495			
0	0.025645146	0.022414174	0.025611562];			
						
						
sigB_sds_inc2 = [0	mean(sigB0) mean(sigB0) mean(sigB0)	
0.01	0.005067031	0.006990135	0.006762378			
0.015	0.012280634	0.009317598	0.015582655			
0.02	0.027419017	0.042910159	0.05185188			
0.025	0.126706769	0.162810843	0.133957153			
0.03	0.153855923	0.195462021	0.190877581];			
						
						
sigB_sds_dec2 = [0.03	0.145835	0.154179935	0.121335801			
0.025	0.122485454	0.12502329	0.124855146			
0.02	0.098273194	0.089598415	0.078894598			
0.015	0.011593509	0.012912144	0.01434043			
0.01	0.008329186	0.009132986	0.009560086			
0	0.008489214	0.009726216	0.009021761];			
						
mprA_sds_inc2 = [0	mean(mprA0) mean(mprA0) mean(mprA0)	
0.01	0.00035919	0.000330751	0.000308034			
0.015	0.000519317	0.000372771	0.000491912			
0.02	0.001581437	0.002607113	0.002256245			
0.025	0.004560364	0.004674801	0.004381904			
0.03	0.006175092	0.006496226	0.006628581];			
						
mprA_sds_dec2 = [0.03	0.007718196	0.00640195	0.005347491			
0.025	0.007741145	0.006191254	0.006599766			
0.02	0.006666336	0.004783149	0.004857503			
0.015	0.000904441	0.000925075	0.000761002			
0.01	0.00077771	0.000645435	0.000694402			
0	0.000654306	0.000712861	0.000827769];

%dynamic error
err(1) = norm(Y(idx_u,6)/Y(1,6) - mean(mprA_0_30(:,2:4),2)/mean(mprA0));%/norm(1 - mean(mprA_0_30(:,2:4),2)/mean(mprA0)); %mprA
err(2) = norm(Y(idx_u,7)/Y(1,7) - mean(sigE_0_30(:,2:4),2)/mean(sigE0));%/norm(1 - mean(sigE_0_30(:,2:4),2)/mean(sigE0)); %sigE
% err(3) = norm(Y(idx_u,13)/Y(1,13)- mean(sigB_0_30(:,2:4),2)/mean(sigB0))/norm(1- mean(sigB_0_30(:,2:4),2)/mean(sigB0)); %sigB
% 
err(4) = norm(X2(idx_d,6)/Y(1,6) - mean(mprA_30_0(:,2:4),2)/mean(mprA0));%/norm(1 - mean(mprA_30_0(:,2:4),2)/mean(mprA0)); %mprA
err(5) = norm(X2(idx_d,7)/Y(1,7) - mean(sigE_30_0(:,2:4),2)/mean(sigE0));%/norm(1 - mean(sigE_30_0(:,2:4),2)/mean(sigE0)); %sigE
% err(6) = norm(X2(idx_d,13)/Y(1,13)- mean(sigB_30_0(:,2:4),2)/mean(sigB0))/norm(1- mean(sigB_30_0(:,2:4),2)/mean(sigB0)); %sigB
% dose response error
err(7) = norm(UP(:,6)/UP(1,6) - mean(mprA_sds_inc2(:,2:4),2)/mean(mprA0));%/norm(1-mean(mprA_sds_inc2(:,2:4),2)/mean(mprA0));
err(8) = norm(UP(:,7)/UP(1,7) - mean(sigE_sds_inc3(:,2:4),2)/mean(sigE0));%/norm(1-mean(sigE_sds_inc3(:,2:4),2)/mean(sigE0));
% err(9) = norm(UP(:,13)/UP(1,13) - mean(sigB_sds_inc2(:,2:4),2)/mean(sigB0))/norm(1-mean(sigB_sds_inc2(:,2:4),2)/mean(sigB0));

err(10) = norm(flipud(Down(:,6))/UP(1,6) - mean(mprA_sds_dec2(:,2:4),2)/mean(mprA0));%/norm(1-mean(mprA_sds_dec2(:,2:4),2)/mean(mprA0));
err(11) = norm(flipud(Down(:,7))/UP(1,7) - mean(sigE_sds_dec3(:,2:4),2)/mean(sigE0));%/norm(1-mean(sigE_sds_dec3(:,2:4),2)/mean(sigE0));
% err(12) = norm(flipud(Down(:,13))/UP(1,13) - mean(sigB_sds_dec2(:,2:4),2)/mean(sigB0))/norm(1-mean(sigB_sds_dec2(:,2:4),2)/mean(sigB0));
% w = [0 0 0 0 0 0 1 1 0 1 1 0]; % ss dose response weights
w = [1 1 0 1 1 0 0 0 0 0 0 0]; % time course weights
%% Error SDS-up mutant deltaE. Data for mprA
out = sum((err.*w).^2);
% else
%     out=1e3;
% end
end