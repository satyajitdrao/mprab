% function out = plotA2_DnaK(x)
% script to simulate MprAB-SigE network in response SDS @modelA2_DnaK called. Output of PSOA2_DnaK may be used
% for determining parameters. plots generated for wild-type models and sigE, mprA deletion mutants
clear global  %comment out if in function mode
clear variables % comment out if in function mode
global ktlnB k2 km2 k3 km3 k4 k5 km5 k6
global kpdeg kmdeg ktlnA ktlnE ktlnR kbtpn1 kbtpn2 kbtpn3 f1 f2 f3 K1 K2 K3 kb kd krdeg
global f5b K5b f1b f3b
global ktlnD kbdf kbdb kdsf kdsb dnakoverex

dnakoverex = 0; % for DnaK overexpression simulation
%% define parameters
% post 06222018
% load('pswarm_A2DnaK_06232018.mat'); x = Solution(1,:);
% x = [-4.2430   -2.9200   -1.6539    1.8510   -5.9753   -5.9993   -6.0000    0.4788    0.8719    1.0516   -3.4816   -1.4059    0.3400 1.6554   -1.4702    1.1186   -2.5518   -3.3592    0.8395   -3.1645   -0.2768   -3.2177   -1.1905   -2.0468   -1.2692   -1.4327 -0.0672   -1.2862   -3.5598   -2.1437    2.3707   -1.4429   -1.4954    0.6437    1.0447];
% load('01042019_ontc2hrdr_dnak.mat'); x = Solution(2,:); % 4 good parameter sets at least
load('pswarm08202019_dnak.mat'); x = Solution(1,:);

% load('09202020_A2.mat'); x = Solution(2,:);
% load('03112019_nokrdeg.mat'); x = Solution(3,:); x(18) = -Inf;
% load('03112019_noultrasens.mat'); x = Solution(3,:);
% load('04132019_nokrdeg.mat'); x = Solution(2,:); x(18) = -Inf;
% load('04132019_noultrasens.mat'); x = Solution(2,:);
% load('05172019_KD_1e-2.mat'); x = Solution(2,:); %load('05172019_KD_1e-1.mat'); load('05172019_KD_5e-2.mat'); load('05172019_KD_1e-2.mat'); 
% load('05172019_krdeg_1e-6.mat'); x = Solution(2,:); %load('05172019_krdeg_1e-5.mat');load('05172019_krdeg_5e-6.mat');load('05172019_krdeg_1e-6.mat');
% load('05172019_NoPhosphatase2'); x = Solution(2,:); x(28:30) = -Inf;

kpdeg= 10^x(1);
kmdeg = 10^x(2);
% translation, transcription
ktlnA = 10^x(3); ktlnA0 = ktlnA;
ktlnB = ktlnA/4.5;
ktlnE = ktlnA/5; ktlnE0 = ktlnE;

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
ktlnD = 10^x(13)*((kbtpn1/kmdeg)*ktlnA);
    % MprA-DNA binding and transcription initiation
f1b = f1; %10^x(36); % f1
f3b = f3; %10^x(37); % f3
K5b = 10^x(16); % 30
    % dose response hill parameters
Kfit1 = 10^x(32);
n1 = 10^x(34);
Kfit2 = 10^x(33);
n2 = 10^x(35);
% post-translational
kdsf = 10^x(15);
kdsb = 10^x(17);
krdeg0 = 10^x(18);
kb = 10^x(19);
kd = 10^x(20);
kbdf = 10^x(21);
kbdb = 10^x(22);
k2 = 10^x(23);
km2 = 10^x(24);
k3= 10^x(25);
km3= 10^x(26);
k4= 10^x(27);
k5 = 10^x(28);
km5 = 10^x(29);
k6 = 10^x(30);
factor2 = 10^x(31);
sdsmax = 10^x(14)*((kbtpn1/(kpdeg*kmdeg))*ktlnA);

% factor2 = 0; % for 09122020 reviewer response
% ktlnE = 0;

% krdeg and Stot as a general function of SDS
krdegf = @(SDS) krdeg0 + krdeg0*factor2*(SDS.^n2./(SDS.^n2 + Kfit2^n2));
v1f = @(SDS) sdsmax*(SDS.^n1./(SDS.^n1 + Kfit1^n1));
%% simulations
X0 = 0*ones(1,16);
% pre-stress steady state
X0([15 16]) = 0; krdeg=krdeg0;
[~, X]=ode15s(@modelA2_DnaK, [0 100*3600], X0); % 0 ss

% % % % total protein levels in pre-stress cells
% SigE = sum(X(end,[3 5]),2) % SigE
% MprA = sum(X(end,[1 2 11 12]),2) %MprA
% MprB = sum(X(end,[8 10 11 12 14]),2) %MprB
% DnaK = sum(X(end,[9 12 14 15]),2) %DnaK
% % % % total protein levels in pre-stress cells

% % krdeg and Stot as a general function of SDS
SDS = 0.03; krdeg = krdegf(SDS); v1 = v1f(SDS);
Y0 = X(end,:); Y0(end,[15 16]) = [0 v1];
tend_up = 2; tend_down = 2;
tdif = 300; tspan = 0:tdif:tend_up*3600; t2hr = (tend_up*3600/tdif)+1;
[tY, Y]=ode15s(@modelA2_DnaK, tspan, Y0); % Post-stress (0 --> 0.03)

X20 = Y(t2hr,:);
X20(9) = Y(t2hr,9) + Y(t2hr,15); % washing step
X20([15 16]) = [0 0]; % washing step
krdeg=krdeg0;
[tX2, X2]=ode15s(@modelA2_DnaK, [0 tend_down]*3600, X20); % 0.03 --> 0
%% dose response simulation
SDS = linspace(0,0.03,100); SDS(2)= []; % range of SDS
UP = zeros(length(SDS),length(X0)); Down = zeros(length(SDS),length(X0));
for i = 1:length(SDS)
% % % % krdeg and Stot as a general function of SDS
Y0(16) = v1f(SDS(i)); X20(16) = v1f(SDS(i));
krdeg = krdegf(SDS(i));
[tU,U] = ode15s(@modelA2_DnaK, [0 tend_up*3600], Y0);
[tD,D] = ode15s(@modelA2_DnaK, [0 tend_down*3600], X20);
UP(i,:) = U(end,:);
Down(i,:) = D(end,:);
end
% 2hr state 0.03->0 and initial state: X(end,:)
UP(1,:) = X(end,:);
Down(1,:) = X2(end,:);
Dtdr = [sum(UP(:,[9 12 14]),2) sum(Down(:,[9 12 14]),2)]; % DT - [D.S]
Btdr = [sum(UP(:,[8 10 11 12 14]),2) sum(Down(:,[8 10 11 12 14]),2)]; %MprB

out.tY = tY; out.Y = Y; out.UP = UP; out.Down = Down;
%% simulate deletion/ disruption mutants
% % delta mprAB
% ktlnA = 0; X0A = 0*ones(1,16);
% X0A([1 2 8 10 11 12 14]) = 0; X0A([15 16]) = 0; krdeg=krdegf(0);
% [~, XA]=ode15s(@modelA2_DnaK, [0 120*3600], X0A); % 0 ss
% % bacteriostatic stress - 2 hrs @ 0.03%
% Y0A = XA(end,:); Y0A(end,[15 16]) = [0 v1f(0.03)]; krdeg = krdegf(0.03);
% [tYA,YA] = ode15s(@modelA2_DnaK, [0 2*3600], Y0A); % 0 --> 0.03
% 
% X2A0 = YA(end,:); 
% X2A0(end,9) = YA(end,9) + YA(end, 15); X2A0(end,[15 16]) = 0; krdeg=krdegf(0);
% [tX2A, X2A]=ode15s(@modelA2_DnaK, [0 2*3600], X2A0); % 0.03 --> 0
% ktlnA = ktlnA0;
% 
% % delta sigE
% X0E = 0*ones(1,16);
% ktlnE = 0; X0E([3 5]) = [0 0];
% X0E([15 16]) = 0; krdeg=krdegf(0);
% [~, XE]=ode15s(@modelA2_DnaK, [0 10*3600], X0E); % 0 ss
% % bacteriostatic stress - 2 hrs @ 0.03%
% krdeg = krdegf(0.03); YE0 = XE(end,:); YE0(end,[15 16]) = [0 v1f(0.03)];
% [tYE,YE] = ode15s(@modelA2_DnaK, [0 2*3600], YE0); % 0 --> 0.03
% X2E0 = YE(end,:); 
% X2E0(end,9) = YE(end,9) + YE(end, 15); X2E0(end,[15 16]) = 0; krdeg=krdegf(0);
% [tX2E, X2E]=ode15s(@modelA2_DnaK, [0 2*3600], X2E0); % 0.03 --> 0
% ktlnE = ktlnE0;
%% graphs
run Rawdata_mprab_sigE_repeat2.m
figure(3);
subplot(2,3,4)
    plot(SDS, UP(:,6)/UP(1,6),'k-'); hold on
    plot(SDS, Down(:,6)/UP(1,6),'r-');
    errorbar(mprA_sds_inc2(:,1),mean(mprA_sds_inc2(:,2:4),2)/mean(mprA0),std(mprA_sds_inc2(:,2:4),'',2)/mean(mprA0),'k^'); hold on
    errorbar(mprA_sds_dec2(:,1),mean(mprA_sds_dec2(:,2:4),2)/mean(mprA0),std(mprA_sds_dec2(:,2:4),'',2)/mean(mprA0),'r^')
    xlabel('SDS (%)'); ylabel('{\it mprA} mRNA (fold)')
    xlim([0 0.035]); ylim([0 15])
subplot(2,3,5)
    plot(SDS, UP(:,7)/UP(1,7),'k-'); hold on
    plot(SDS, Down(:,7)/UP(1,7),'r-');
    errorbar(sigE_sds_inc3(:,1),mean(sigE_sds_inc3(:,2:4),2)/mean(sigE0),std(sigE_sds_inc3(:,2:4),'',2)/mean(sigE0),'k^'); hold on
    errorbar(sigE_sds_dec3(:,1),mean(sigE_sds_dec3(:,2:4),2)/mean(sigE0),std(sigE_sds_dec3(:,2:4),'',2)/mean(sigE0),'r^')
    xlim([0 0.035]); ylim([0 15])       
    xlabel('SDS (%)'); ylabel('{\it sigE} mRNA (fold)')
subplot(2,3,6)
    plot(SDS, UP(:,13)/UP(1,13),'k-'); hold on
    plot(SDS, Down(:,13)/UP(1,13),'r-');
    errorbar(sigB_sds_inc2(:,1),mean(sigB_sds_inc2(:,2:4),2)/mean(sigB0),std(sigB_sds_inc2(:,2:4),'',2)/mean(sigB0),'k^'); hold on
    errorbar(sigB_sds_dec2(:,1),mean(sigB_sds_dec2(:,2:4),2)/mean(sigB0),std(sigB_sds_dec2(:,2:4),'',2)/mean(sigB0),'r^')
    xlabel('SDS'); 
    xlim([0 0.035]); ylim([0 30])
subplot(2,3,1)
    errorbar(mprA_0_30(:,1),mean(mprA_0_30(:,2:4),2)/mean(mprA0),std(mprA_0_30(:,2:4),'',2)/mean(mprA0),'k^'); hold on
    errorbar(mprA_30_0(:,1)+120,mean(mprA_30_0(:,2:4),2)/mean(mprA0),std(mprA_30_0(:,2:4),'',2)/mean(mprA0),'r^')
    plot(tY/60, Y(:,6)/X(end,6),'k', tX2/60 +120, X2(:,6)/X(end,6),'r'); hold on
%     plot(tYA/60, YA(:,6)/X(end,6),'k:', tX2A/60 +120, X2A(:,6)/X(end,6),'r:');
%     plot(tYE/60, YE(:,6)/X(end,6),'k-.', tX2E/60 +120, X2E(:,6)/X(end,6),'r-.');
    xlabel('Time (mins)'); ylabel('{\it mprA} mRNA (fold)')
    xlim([-10 250]); ylim([0 12])
subplot(2,3,2)
    errorbar(sigE_0_30(:,1),mean(sigE_0_30(:,2:4),2)/mean(sigE0),std(sigE_0_30(:,2:4),'',2)/mean(sigE0),'k^'); hold on
    errorbar(sigE_30_0(:,1)+120,mean(sigE_30_0(:,2:4),2)/mean(sigE0),std(sigE_30_0(:,2:4),'',2)/mean(sigE0),'r^')
    plot(tY/60, Y(:,7)/X(end,7),'k', tX2/60 +120, X2(:,7)/X(end,7),'r'); hold on
%     plot(tYA/60, YA(:,7)/X(end,7),'k:', tX2A/60 +120, X2A(:,7)/X(end,7),'r:');
%     plot(tYE/60, YE(:,7)/X(end,7),'k-.', tX2E/60 +120, X2E(:,7)/X(end,7),'r-.');
    xlabel('Time (mins)'); ylabel('{\it sigE} mRNA (fold)')
    xlim([-10 250]); ylim([0 15])
subplot(2,3,3)
    errorbar(sigB_0_30(:,1),mean(sigB_0_30(:,2:4),2)/mean(sigB0),std(sigB_0_30(:,2:4),'',2)/mean(sigB0),'k^'); hold on
    errorbar(sigB_30_0(:,1)+120,mean(sigB_30_0(:,2:4),2)/mean(sigB0),std(sigB_30_0(:,2:4),'',2)/mean(sigB0),'r^')
    plot(tY/60, Y(:,13)/X(end,13),'k', tX2/60 +120, X2(:,13)/X(end,13),'r'); hold on
%     plot(tYA/60, YA(:,13)/X(end,13),'k:', tX2A/60 +120, X2A(:,13)/X(end,13),'r:');
%     plot(tYE/60, YE(:,13)/X(end,13),'k-.', tX2E/60 +120, X2E(:,13)/X(end,13),'r-.');
    xlabel('time(min)')
    title('sigB','FontSize',16)
    xlim([-10 250]); ylim([0 30])
set(findobj(gcf,'Type','axes'), 'linewidth',1,'fontname','arial','fontsize',15)    
%% parameter set distributions
% pars = 'kpdeg kmdeg ktlnA ktlnB ktlnE ktlnR kbtpn1 kbtpn2 kbtpn3 f1   f2 f3  K1   K2   Kfit2  ktlnD SDS   kdsf K5b  Kfit1  kdsb n1  krdeg kb kd   kbdf kbdb k2   km2 k3 km3  k4  k5   km5  k6    factor2 n2';
% pars = 'kpdeg kmdeg ktlnA ktlnR kbtpn1 kbtpn2 kbtpn3 f1   f2  f3  K1   K2   ktlnD SDS  kdsf K5b  kdsb krdeg kb kd   kbdf kbdb k2   km2 k3 km3  k4  k5   km5  k6    factor2  kfit1 kfit2 n1  n2';
% parindex = strsplit(pars,' ');
% sol = Solution(:,:); x = Solution(1,:);
% figure%(2)
% for i = 1:size(sol,2)-1
%     if sol(end-1,i) ~= sol(end,i)
% subplot(6,6,i)
%     hold on
%     [n, y] = hist(sol(1:end-2,i));
%     h = bar(y,n,'hist');
%     set(h,'FaceColor',[1 0 1],'facealpha',.5,'edgecolor',[1 1 1])
%     set(gca,'XLim',sol([end-1 end],i))
% plot(x(i)*ones(1,20),linspace(0,10,20),'linewidth',1.8,'color','k')
%     xlabel(parindex{i})
%     end
% end

%% dnak overexpression: pmprA-dnak strain
dnakoverex = 1; % for DnaK overexpression simulation
% simulations
X0 = 0*ones(1,16); ktlnE=0;
% pre-stress steady state
X0([15 16]) = 0; krdeg=krdeg0;
[~, X]=ode15s(@modelA2_DnaK, [0 100*3600], X0); % 0 ss
SDS = 0.03; krdeg = krdegf(SDS); v1 = v1f(SDS);
tend_up = 2; tend_down = 2;
tdif = 300; tspan = 0:tdif:tend_up*3600; t2hr = (tend_up*3600/tdif)+1;
Y0 = X(end,:); Y0(end,[15 16]) = [0 v1];
[tY, Y]=ode15s(@modelA2_DnaK, tspan, Y0); % Post-stress (0 --> 0.03)

X20 = Y(t2hr,:);
X20(9) = Y(t2hr,9) + Y(t2hr,15); % washing step
X20([15 16]) = [0 0]; % washing step
krdeg=krdeg0;
[tX2, X2]=ode15s(@modelA2_DnaK, [0 tend_down]*3600, X20); % 0.03 --> 0
% dose response simulation
SDS = linspace(0,0.03,100); SDS(2)= []; % range of SDS
UP = zeros(length(SDS),length(X0)); Down = zeros(length(SDS),length(X0));
for i = 1:length(SDS)
% % % % krdeg and Stot as a general function of SDS
Y0(16) = v1f(SDS(i)); X20(16) = v1f(SDS(i));
krdeg = krdegf(SDS(i));
[tU,U] = ode15s(@modelA2_DnaK, [0 tend_up*3600], Y0);
[tD,D] = ode15s(@modelA2_DnaK, [0 tend_down*3600], X20);
UP(i,:) = U(end,:);
Down(i,:) = D(end,:);
end
% 2hr state 0.03->0 and initial state: X(end,:)
UP(1,:) = X(end,:);
Down(1,:) = X2(end,:);
Dtdr = [sum(UP(:,[9 12 14]),2) sum(Down(:,[9 12 14]),2)]; % DT - [D.S]
Btdr = [sum(UP(:,[8 10 11 12 14]),2) sum(Down(:,[8 10 11 12 14]),2)]; %MprB

out.tY = tY; out.Y = Y; out.UP = UP; out.Down = Down;

sigeup = [0	0.009222857	0.006911959	0.006243846
0.01	0.003868277	0.004513294	0.006778879
0.015	0.003227023	0.006806088	0.006757394
0.02	0.006760242	0.009335605	0.004931927
0.025	0.006811719	0.008471673	0.008489242
0.03	0.010218423	0.010344646	0.01000642];
			
sigedown=[0.03	0.010218423	0.010344646	0.01000642
0.025	0.010723788	0.008858876	0.00850916
0.02	0.006475443	0.007851588	0.012831686
0.015	0.004340508	0.008212714	0.007799776
0.01	0.004262206	0.004598128	0.00358271
0	0.003706889	0.005501234	0.005767781];

sigbup = [0	0.00409905	0.004156425	0.00542109
0.01	0.004069337	0.003845107	0.002987468
0.015	0.002767439	0.004006565	0.006825838
0.02	0.003742056	0.005293503	0.00397242
0.025	0.007451907	0.008383343	0.007505207
0.03	0.007681679	0.008419852	0.008288809];
			
sigbdown =[0.03	0.007681679	0.008419852	0.008288809
0.025	0.006303314	0.006694711	0.006077797
0.02	0.006694711	0.006470452	0.007578441
0.015	0.006077797	0.006446056	0.008014051
0.01	0.002229145	0.003045959	0.002444401
0	0.002447361	0.002132161	0.003469446];
			
mpraup=[0	0.001514878	0.001091251	0.001853778
0.01	0.001212741	0.001248998	0.00128021
0.015	0.00093946	0.001184002	0.001432977
0.02	0.001715592	0.002490786	0.001393028
0.025	0.003650811	0.003869666	0.003011294
0.03	0.004366534	0.004167552	0.003568262];
			
mpradown =[0.03	0.004366534	0.004167552	0.003568262
0.025	0.004286337	0.003908929	0.003157544
0.02	0.002131914	0.00244504	0.00320252
0.015	0.001297892	0.002518246	0.001217401
0.01	0.000687199	0.000892745	0.000764871
0	0.000551388	0.000888415	0.000813804];

mprA0 = mpraup(1,(2:4));
sigE0= sigeup(1,(2:4));
sigB0= sigbup(1,(2:4));

figure(10);
subplot(1,2,1)
    errorbar(mpraup(:,1),mean(mpraup(:,2:4),2)/mean(mprA0),std(mpraup(:,2:4),'',2)/mean(mprA0),'k^', 'linewidth',2); hold on
    errorbar(mpradown(:,1),mean(mpradown(:,2:4),2)/mean(mprA0),std(mpradown(:,2:4),'',2)/mean(mprA0),'r^', 'linewidth',2)
    plot(SDS, UP(:,6)/UP(1,6),'k-', 'linewidth',2); hold on
    plot(SDS, Down(:,6)/UP(1,6),'r-', 'linewidth',2);
    xlabel('SDS'); ylabel('{\it mprA} mRNA (fold)');
%     xlim([0 0.035]); ylim([0 15])
subplot(1,2,2)
    errorbar(sigeup(:,1),mean(sigeup(:,2:4),2)/mean(sigE0),std(sigeup(:,2:4),'',2)/mean(sigE0),'k^', 'linewidth',2); hold on
    errorbar(sigedown(:,1),mean(sigedown(:,2:4),2)/mean(sigE0),std(sigedown(:,2:4),'',2)/mean(sigE0),'r^', 'linewidth',2)
    plot(SDS, UP(:,7)/UP(1,7),'k-', 'linewidth',2); hold on
    plot(SDS, Down(:,7)/UP(1,7),'r-', 'linewidth',2);
    xlabel('SDS'); ylabel('{\it sigE} mRNA (fold)');
set(findobj(gcf,'Type','axes'), 'linewidth',1.5,'fontname','arial','fontsize',15,'XLim',[0 0.03],'YLim',[0 5])

    
figure(11) %09122020 msystems reviewer #2 point 10 analysis 
subplot(2,3,1);
plot(tY/3600, Y(:,2),'linewidth',2); hold on;
xlabel('time (hrs)'); ylabel('MprA-P (\muM)'); %ylim([0 7e-3])
subplot(2,3,2);
plot(tY/3600, Y(:,1), 'linewidth',2); hold on;
xlabel('time (hrs)'); ylabel('MprA (\muM)'); %ylim([0 6])
subplot(2,3,3);
plot(tY/3600, Y(:,6)/Y(1,6), 'linewidth',2); hold on;
xlabel('time (hrs)'); ylabel('{\it mprA} mRNA (fold)'); %ylim([0 7])
subplot(2,3,4);
plot(tY/3600, Y(:,7)/Y(1,7), 'linewidth',2); hold on;
xlabel('time (hrs)'); ylabel('{\it sigE} mRNA (fold)'); %ylim([0 7])
subplot(2,3,5);
plot(tY/3600, sum(Y(:,[1 2 11 12]),2)/sum(X(end,[1 2 11 12]),2), 'linewidth',2); hold on;
xlabel('time (hrs)'); ylabel('Total-MprA (fold)'); %ylim([0 7])
set(findobj(gcf,'Type','axes'), 'linewidth',1.5,'fontname','arial','fontsize',15,'fontweight','bold')
dnakoverex = 0; % for DnaK overexpression simulation