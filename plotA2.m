function out = plotA2(x) % comment out in default mode
% simulations of modelA2 (without DnaK) using parameter sets
% from PSOA2. Time-course and dose-response simulations plotted against
% experimental data for mprA, sigE.
% clear global % comment out in function mode
% clear variables % comment out in function mode
global ktlnB k1 km1 k2 km2 k3 km3 k4 k5 km5 k6
global kpdeg kmdeg ktlnA ktlnE ktlnR kbtpn1 kbtpn2 kbtpn3 f1 f2 f3 K1 K2 K3 kb kd krdeg
global f5b K5b

% load('06222018_modelA2_dr_1.mat') % only 1 parameter set, fits dose response (8 hr)
% load('pswarm_A2_dr_06222018.mat'); x = Solution(2,:); % only 2 fits dr (8 hr)
% load('12192018_ontc_ssdr.mat'); x = Solution(3,:); % fits ON time course only

% important for paper:
load('06232018_modelA2_dr_1') % 1 parameter set, fits true steady state dose response w/o A2(used in paper figure, fmincon used to get this)
% load('12242018_ssdr_nodnak.mat'); x = Solution(3,:); % fits dose response (true steady state, w/o A2)
% load('modelA2_0120_rescaledlimits'); x = Solution(1,:); % fits all time courses (w/A2)
% load('08222019_modelA2_tconly.mat'); x = Solution(25,:);
% load('08222019_modelA2_DRonly.mat'); x = Solution(4,:);

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
K5b =  Inf; %10^x(13);%
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
%% simulations
krdegf = @(SDS) krdeg0 + krdeg0*factor2*(SDS.^n2./(SDS.^n2 + Kfit2^n2));
k1f = @(SDS) k10+k10*factor1*(SDS.^n1./(SDS.^n1 + Kfit1^n1));
%      A  Ap   E   R  ER  mAB   mE    Bs  B   Bp  kpa   pap   sigB
X0 = zeros(1,13); %[0.1 0.1 0.1 0.1 0.1 0.001 0.001 0.1 0.1 0.1 0.001 0.001 0.001];
tspan = [0 2]*3600;
% tspan = 0:5*3600:180*3600;
% pre-stress steady state
k1=k10;
krdeg=krdeg0;
[tX, X]=ode15s(@modelA2, [0 100*3600], X0); % 0 ss
% bacteriostatic stress - 2 hrs @ 0.03%
k1 = k1f(0.03);% k1 = k10*factor1;
krdeg = krdegf(0.03);% krdeg = krdeg0*factor2;
[tY,Y] = ode15s(@modelA2, tspan, X(end,:)); % 0 --> 0.03
% 0.03 --> 0% SDS
k1=k10;
krdeg=krdeg0;
[tX2, X2]=ode15s(@modelA2, tspan, Y(end,:)); % 0.03 --> 0

SDS = linspace(0,0.03,37); SDS(2) = [];

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

out.Y = Y; out.UP = UP; out.Down = Down; out.tY = tY; % comment out in default mode
%% simulate deletion/ disruption mutants
% % % delta mprAB
% % ktlnA = 0; X0([1 2 8:12]) = zeros(1,length([1 2 8:12]));
% % k1=k10;
% % krdeg=krdeg0;
% % [~, XA]=ode15s(@modelA2, [0 10*3600], X0); % 0 ss
% % % bacteriostatic stress - 2 hrs @ 0.03%
% % k1 = k1f(0.03); %k10*factor1;
% % krdeg = krdegf(0.03); %krdeg0*factor2;
% % [tYA,YA] = ode15s(@modelA2, [0 2*3600], XA(end,:)); % 0 --> 0.03
% % k1=k10;
% % krdeg=krdeg0;
% % [tX2A, X2A]=ode15s(@modelA2, [0 2*3600], YA(end,:)); % 0.03 --> 0
% % X0 = 0.01*ones(1,13);
% % ktlnA = ktlnA0;
% % 
% % % delta sigE
% % X0 = 0.01*ones(1,13);
% % ktlnE = 0; X0([3 5]) = [0 0];
% % k1=k10;
% % krdeg=krdeg0;
% % [~, XE]=ode15s(@modelA2, [0 10*3600], X0); % 0 ss
% % % bacteriostatic stress - 2 hrs @ 0.03%
% % k1 = k1f(0.03); % k10*factor1;
% % krdeg = krdegf(0.03); %krdeg0*factor2;
% % [tYE,YE] = ode15s(@modelA2, [0 2*3600], XE(end,:)); % 0 --> 0.03
% % k1=k10;
% % krdeg=krdeg0;
% % [tX2E, X2E]=ode15s(@modelA2, [0 2*3600], YE(end,:)); % 0.03 --> 0
% % ktlnE = ktlnE0;
% % X0 = 0.01*ones(1,13);
%% plotting
% % run Rawdata_mprab_sigE_repeat2.m
% % figure(1);
% % subplot(2,3,4)
% % plot(SDS, UP(:,6)/UP(1,6),'k'); hold on
% % plot(SDS, Down(:,6)/UP(1,6),'r');
% % errorbar(mprA_sds_inc2(:,1),mean(mprA_sds_inc2(:,2:4),2)/mean(mprA0),std(mprA_sds_inc2(:,2:4),'',2)/mean(mprA0),'k^'); hold on
% % errorbar(mprA_sds_dec2(:,1),mean(mprA_sds_dec2(:,2:4),2)/mean(mprA0),std(mprA_sds_dec2(:,2:4),'',2)/mean(mprA0),'r^')
% % axis([0 0.035 0 15])
% % xlabel('SDS (%)');
% % ylabel('{\it mprA}(fold-change)');
% % legend('Increase','Decrease','location','northwest'); legend boxoff
% % subplot(2,3,5)
% % plot(SDS, UP(:,7)/UP(1,7),'k'); hold on
% % plot(SDS, Down(:,7)/UP(1,7),'r');
% % errorbar(sigE_sds_inc3(:,1),mean(sigE_sds_inc3(:,2:4),2)/mean(sigE0),std(sigE_sds_inc3(:,2:4),'',2)/mean(sigE0),'k^'); hold on
% % errorbar(sigE_sds_dec3(:,1),mean(sigE_sds_dec3(:,2:4),2)/mean(sigE0),std(sigE_sds_dec3(:,2:4),'',2)/mean(sigE0),'r^')
% % axis([0 0.035 0 15])
% % ylabel('{\it sigE}(fold-change)');
% % xlabel('SDS (%)');
% % subplot(2,3,6)
% % plot(SDS, UP(:,13)/UP(1,13),'k'); hold on
% % plot(SDS, Down(:,13)/UP(1,13),'r');
% % errorbar(sigB_sds_inc2(:,1),mean(sigB_sds_inc2(:,2:4),2)/mean(sigB0),std(sigB_sds_inc2(:,2:4),'',2)/mean(sigB0),'k^'); hold on
% % errorbar(sigB_sds_dec2(:,1),mean(sigB_sds_dec2(:,2:4),2)/mean(sigB0),std(sigB_sds_dec2(:,2:4),'',2)/mean(sigB0),'r^')
% % axis([0 0.035 0 26])
% % ylabel('{\it sigB}(fold-change)'); 
% % xlabel('SDS (%)');
% % figure(1);
% %     subplot(2,3,1)
% %         errorbar(mprA_0_30(:,1),mean(mprA_0_30(:,2:4),2)/mean(mprA0),std(mprA_0_30(:,2:4),'',2)/mean(mprA0),'k^'); hold on
% %         errorbar(mprA_30_0(:,1)+120,mean(mprA_30_0(:,2:4),2)/mean(mprA0),std(mprA_30_0(:,2:4),'',2)/mean(mprA0),'r^')
% %         plot(tY/60, Y(:,6)/X(end,6),'k', tX2/60 +120, X2(:,6)/X(end,6),'r'); hold on
% % %         plot(tYA/60, YA(:,6)/X(end,6),'k:', tX2A/60 +120, X2A(:,6)/X(end,6),'r:');
% % %         plot(tYE/60, YE(:,6)/X(end,6),'k-.', tX2E/60 +120, X2E(:,6)/X(end,6),'r-.');
% % ylabel('{\it mprA}(fold-change)'); xlabel('time (mins)');
% %         axis([-10 250 0 15])
% %         legend('SDS 0 \rightarrow 0.03%','SDS 0.03 \rightarrow 0%','location','northwest'); legend boxoff
% %     subplot(2,3,2)
% %         errorbar(sigE_0_30(:,1),mean(sigE_0_30(:,2:4),2)/mean(sigE0),std(sigE_0_30(:,2:4),'',2)/mean(sigE0),'k^'); hold on
% %         errorbar(sigE_30_0(:,1)+120,mean(sigE_30_0(:,2:4),2)/mean(sigE0),std(sigE_30_0(:,2:4),'',2)/mean(sigE0),'r^')
% %         plot(tY/60, Y(:,7)/X(end,7),'k', tX2/60 +120, X2(:,7)/X(end,7),'r'); hold on
% % %         plot(tYA/60, YA(:,7)/X(end,7),'k:', tX2A/60 +120, X2A(:,7)/X(end,7),'r:');
% % %         plot(tYE/60, YE(:,7)/X(end,7),'k-.', tX2E/60 +120, X2E(:,7)/X(end,7),'r-.');
% % ylabel('{\it sigE}(fold-change)'); xlabel('time (mins)');
% %         axis([-10 250 0 15])
% %     subplot(2,3,3)
% %         errorbar(sigB_0_30(:,1),mean(sigB_0_30(:,2:4),2)/mean(sigB0),std(sigB_0_30(:,2:4),'',2)/mean(sigB0),'k^'); hold on
% %         errorbar(sigB_30_0(:,1)+120,mean(sigB_30_0(:,2:4),2)/mean(sigB0),std(sigB_30_0(:,2:4),'',2)/mean(sigB0),'r^')
% %         plot(tY/60, Y(:,13)/X(end,13),'k', tX2/60 +120, X2(:,13)/X(end,13),'r'); hold on
% % %         plot(tYA/60, YA(:,13)/X(end,13),'k:', tX2A/60 +120, X2A(:,13)/X(end,13),'r:');
% % %         plot(tYE/60, YE(:,13)/X(end,13),'k-.', tX2E/60 +120, X2E(:,13)/X(end,13),'r-.');
% % ylabel('{\it sigB}(fold-change)'); xlabel('time (mins)');
% %         axis([-10 250 0 26])
% % set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',13, ...
% % 'LineWidth', 1,'layer','top');
% % set(gcf,'position',[10,10,468,250]) 


% % param = 'kpdeg kmdeg ktlnA ktlnB ktlnE ktlnR kbtpn1 kbtpn2 kbtpn3 f1 f2  f3  K1 K2 K3 f5a f5b K5a K5b f6 K4 kb4 krdeg kb kd k1 km1 k2 km2 k3 km3 k4 k5 km5 k6 factor1 factor2';
% % parindex = strsplit(param, ' ');
% % 
% % figure(2)
% % for i = 1:size(sol,2)-2
% %     if sol(end-1,i) ~= sol(end,i)
% %     subplot(5,8,i)
% %     hold on
% %     [n, x] = hist(sol(1:end-2,i));
% %     h = bar(x,n,'hist');
% %     set(h,'FaceColor',[0 0 1],'facealpha',.5,'edgecolor',[1 1 1])
% %     set(gca,'XLim',sol([end-1 end],i))
% %     
% %     xlabel(parindex{i})
% %     end
% % end


% fmincon(@errA2,x,[],[],[],[],lb,ub,[],options2)


%% paper figure
% % % parameter type I
% % run Rawdata_mprab_sigE_repeat2.m
% % 
% % figure(1);
% % subplot(2,4,5)
% % plot(SDS, UP(:,6)/UP(1,6),'k'); hold on
% % plot(SDS, Down(:,6)/UP(1,6),'r');
% % errorbar(mprA_sds_inc2(:,1),mean(mprA_sds_inc2(:,2:4),2)/mean(mprA0),std(mprA_sds_inc2(:,2:4),'',2)/mean(mprA0),'k^'); hold on
% % errorbar(mprA_sds_dec2(:,1),mean(mprA_sds_dec2(:,2:4),2)/mean(mprA0),std(mprA_sds_dec2(:,2:4),'',2)/mean(mprA0),'r^')
% % axis([0 0.035 0 15])
% % xlabel('SDS (%)');
% % ylabel('{\it mprA}(fold-change)');
% % legend('Increase','Decrease','location','northwest'); legend boxoff
% % subplot(2,4,6)
% % plot(SDS, UP(:,7)/UP(1,7),'k'); hold on
% % plot(SDS, Down(:,7)/UP(1,7),'r');
% % errorbar(sigE_sds_inc3(:,1),mean(sigE_sds_inc3(:,2:4),2)/mean(sigE0),std(sigE_sds_inc3(:,2:4),'',2)/mean(sigE0),'k^'); hold on
% % errorbar(sigE_sds_dec3(:,1),mean(sigE_sds_dec3(:,2:4),2)/mean(sigE0),std(sigE_sds_dec3(:,2:4),'',2)/mean(sigE0),'r^')
% % axis([0 0.035 0 15])
% % ylabel('{\it sigE}(fold-change)');
% % xlabel('SDS (%)');
% %         subplot(2,4,1)
% %         errorbar(mprA_0_30(:,1),mean(mprA_0_30(:,2:4),2)/mean(mprA0),std(mprA_0_30(:,2:4),'',2)/mean(mprA0),'k^'); hold on
% %         errorbar(mprA_30_0(:,1)+120,mean(mprA_30_0(:,2:4),2)/mean(mprA0),std(mprA_30_0(:,2:4),'',2)/mean(mprA0),'r^')
% %         plot(tY/60, Y(:,6)/X(end,6),'k', tX2/60 +120, X2(:,6)/X(end,6),'r'); hold on
% %         ylabel('{\it mprA}(fold-change)'); xlabel('time (mins)');
% %         axis([-10 250 0 15])
% %         legend('SDS 0 \rightarrow 0.03%','SDS 0.03 \rightarrow 0%','location','northwest'); legend boxoff
% %         subplot(2,4,2)
% %         errorbar(sigE_0_30(:,1),mean(sigE_0_30(:,2:4),2)/mean(sigE0),std(sigE_0_30(:,2:4),'',2)/mean(sigE0),'k^'); hold on
% %         errorbar(sigE_30_0(:,1)+120,mean(sigE_30_0(:,2:4),2)/mean(sigE0),std(sigE_30_0(:,2:4),'',2)/mean(sigE0),'r^')
% %         plot(tY/60, Y(:,7)/X(end,7),'k', tX2/60 +120, X2(:,7)/X(end,7),'r'); hold on
% %         ylabel('{\it sigE}(fold-change)'); xlabel('time (mins)');
% %         axis([-10 250 0 15])
% % % parameter type II
% % figure(1);
% % subplot(2,4,7)
% % plot(SDS, UP(:,6)/UP(1,6),'k'); hold on
% % plot(SDS, Down(:,6)/UP(1,6),'r');
% % errorbar(mprA_sds_inc2(:,1),mean(mprA_sds_inc2(:,2:4),2)/mean(mprA0),std(mprA_sds_inc2(:,2:4),'',2)/mean(mprA0),'k^'); hold on
% % errorbar(mprA_sds_dec2(:,1),mean(mprA_sds_dec2(:,2:4),2)/mean(mprA0),std(mprA_sds_dec2(:,2:4),'',2)/mean(mprA0),'r^')
% % axis([0 0.035 0 15])
% % xlabel('SDS (%)');
% % ylabel('{\it mprA}(fold-change)');
% % legend('Increase','Decrease','location','northwest'); legend boxoff
% % subplot(2,4,8)
% % plot(SDS, UP(:,7)/UP(1,7),'k'); hold on
% % plot(SDS, Down(:,7)/UP(1,7),'r');
% % errorbar(sigE_sds_inc3(:,1),mean(sigE_sds_inc3(:,2:4),2)/mean(sigE0),std(sigE_sds_inc3(:,2:4),'',2)/mean(sigE0),'k^'); hold on
% % errorbar(sigE_sds_dec3(:,1),mean(sigE_sds_dec3(:,2:4),2)/mean(sigE0),std(sigE_sds_dec3(:,2:4),'',2)/mean(sigE0),'r^')
% % axis([0 0.035 0 15])
% % ylabel('{\it sigE}(fold-change)');
% % xlabel('SDS (%)');
% %         subplot(2,4,3)
% %         errorbar(mprA_0_30(:,1),mean(mprA_0_30(:,2:4),2)/mean(mprA0),std(mprA_0_30(:,2:4),'',2)/mean(mprA0),'k^'); hold on
% %         errorbar(mprA_30_0(:,1)+120,mean(mprA_30_0(:,2:4),2)/mean(mprA0),std(mprA_30_0(:,2:4),'',2)/mean(mprA0),'r^')
% %         plot(tY/60, Y(:,6)/X(end,6),'k', tX2/60 +120, X2(:,6)/X(end,6),'r'); hold on
% %         ylabel('{\it mprA}(fold-change)'); xlabel('time (mins)');
% %         axis([-10 250 0 15])
% %         legend('SDS 0 \rightarrow 0.03%','SDS 0.03 \rightarrow 0%','location','northwest'); legend boxoff
% %         subplot(2,4,4)
% %         errorbar(sigE_0_30(:,1),mean(sigE_0_30(:,2:4),2)/mean(sigE0),std(sigE_0_30(:,2:4),'',2)/mean(sigE0),'k^'); hold on
% %         errorbar(sigE_30_0(:,1)+120,mean(sigE_30_0(:,2:4),2)/mean(sigE0),std(sigE_30_0(:,2:4),'',2)/mean(sigE0),'r^')
% %         plot(tY/60, Y(:,7)/X(end,7),'k', tX2/60 +120, X2(:,7)/X(end,7),'r'); hold on
% %         ylabel('{\it sigE}(fold-change)'); xlabel('time (mins)');
% %         axis([-10 250 0 15])
% %         
% %         set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',13, ...
% % 'LineWidth', 1,'layer','top');
% % set(gcf,'position',[10,10,468,250]) 