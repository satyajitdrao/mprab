% % % Experiments
run Rawdata_mprab_sigE_repeat2.m
expmpraup = mean(mprA_sds_inc2(:,2:4),2)/mean(mprA0);
expmpradn = flipud(mean(mprA_sds_dec2(:,2:4),2)/mean(mprA0));
stdexpmpraup = std(mprA_sds_inc2(:,2:4),'',2)/mean(mprA0);
stdexpmpradn = std(mprA_sds_dec2(:,2:4),'',2)/mean(mprA0);
sqrt((stdexpmpraup.^2)/3+(stdexpmpradn.^2)/3)
expsigeup = mean(sigE_sds_inc3(:,2:4),2)/mean(sigE0);
expsigedn = flipud(mean(sigE_sds_dec3(:,2:4),2)/mean(sigE0));
stdexpsigeup = std(sigE_sds_inc3(:,2:4),'',2)/mean(sigE0);
stdexpsigedn = std(sigE_sds_dec3(:,2:4),'',2)/mean(sigE0);
hexpstd = sqrt((stdexpsigeup.^2)/3+(stdexpsigedn.^2)/3);
% dose-response
% measure for hysteresis: experiments
hexp(1) = mean((expmpradn([4 5]) - expmpraup([4 5]))); %mpra
hexp(2) = mean((expsigedn([4 5]) - expsigeup([4 5]))); %sige
% time course response time
tmpra = mprA_0_30(:,1);
mpraup = mean(mprA_0_30(:,2:4),2);
spline3 = fit(tmpra,mpraup, 'smoothingspline');
tdif = 5;
tau90exp(1) = ((find(spline3(0:tdif:2*60) > spline3(0)+ 0.9*(spline3(2*60)-spline3(0)),1,'first'))-1)*tdif;
tsige = sigE_0_30(:,1);
sigeup = mean(sigE_0_30(:,2:4),2);
spline1 = fit(tsige,sigeup, 'smoothingspline');
tau90exp(2)=((find(spline1(0:tdif:2*60) > spline1(0)+ 0.9*(spline1(2*60)-spline1(0)),1,'first'))-1)*tdif;


% Simulations
load('superset non dnak time course')
% load('superset non dnak hysteresis')
% load('superset dnak')
n = size(supset,1);
for j = 1:n
    x = supset(j,:);
    out = plotA2(x); Down = out.Down; Y = out.Y; UP = out.UP; tY = out.tY;
%     out = plotA2_DnaK(x); Down = out.Down; Y = out.Y; UP = out.UP; tY = out.tY;
% Measure for hysteresis
hsim(j,1) = mean((Down([4 5],6)/UP(1,6) - UP([4 5],6)/UP(1,6)));
hsim(j,2) = mean((Down([4 5],7)/UP(1,7) - UP([4 5],7)/UP(1,7)));
% t90
tau90(j,1) = tY(find(Y(:,6)> Y(1,6)+ 0.9*(Y(end,6)-Y(1,6)),1,'first'))/60;
tau90(j,2) = tY(find(Y(:,7)> Y(1,7)+ 0.9*(Y(end,7)-Y(1,7)),1,'first'))/60;
end
quant(:,:,1) = hsim; quant(:,:,2) = tau90;
% % % goodness of fit - dose-response
% % gof(1) = 1-norm(UP(:,6)/UP(1,6) - mean(mprA_sds_inc2(:,2:4),2)/mean(mprA0))^2/norm(1-mean(mprA_sds_inc2(:,2:4),2)/mean(mprA0))^2;
% % gof(2) = 1-norm(UP(:,7)/UP(1,7) - mean(sigE_sds_inc3(:,2:4),2)/mean(sigE0))^2/norm(1-mean(sigE_sds_inc3(:,2:4),2)/mean(sigE0))^2;
% % gof(3) = 1-norm(flipud(Down(:,6))/UP(1,6) - mean(mprA_sds_dec2(:,2:4),2)/mean(mprA0))^2/norm(1-mean(mprA_sds_dec2(:,2:4),2)/mean(mprA0))^2;
% % gof(4) = 1-norm(flipud(Down(:,7))/UP(1,7) - mean(sigE_sds_dec3(:,2:4),2)/mean(sigE0))^2/norm(1-mean(sigE_sds_dec3(:,2:4),2)/mean(sigE0))^2;


hsim = quant(:,:,1); tau90 = quant(:,:,2);
figure(6); loglog(hsim(:,2),tau90(:,2),'kd')
hold on; loglog(hexp(2)*ones(1,20),logspace(1,4,20),'k')
loglog((hexp(2)+0.5158)*ones(1,20),logspace(1,4,20),'k');loglog((hexp(2)-0.5158)*ones(1,20),logspace(1,4,20),'k')
loglog(logspace(-3,1,20),tau90exp(2)*ones(1,20),'k')
loglog(logspace(-3,1,20),15*ones(1,20),'k');loglog(logspace(-3,1,20),30*ones(1,20),'k');
xlim([1e-3 1e1]); ylim([1e1 1e4])
ylabel('Response time (sigE) \tau90(min)');
xlabel('Hysteresis measure (sigE)')
legend('Parameter type I','Parameter type II','experimental hysteresis','experimental response time')
% figure(7); semilogy(hsim(2),tau90(2),'mo')
% hold on; semilogy(hexp(2),tau90exp(2),'kd')

figure; semilogx([1e0 1e-1 1e-2],[95 85 65],'k--o'); xlabel('K_D'); ylabel('mprA mRNA \tau90(min)'); hold on; semilogx(logspace(-2,0,10), 60*ones(1,10),'k','linewidth',2)
figure; semilogx([1e-4 1e-5 5e-6 1e-6],[20 55 65 95],'k--o'); xlabel('RseA degradation rate (s^{-1}'); ylabel('mprA mRNA \tau50(min)'); hold on; semilogx(logspace(-6,-4,10), 20*ones(1,10),'k','linewidth',2)