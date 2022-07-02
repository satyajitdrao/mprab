run Rawdataplots.m
run Rawdata_050417.m
run Rawdata_mprab_sigE_repeat2.m
figure(1);
    subplot(2,3,1)
        errorbar(mprA_wt(:,1),mprA_wt(:,2)/mprA_wt(1,2),mprA_wt_std/mprA_wt(1,2),'k-*')
        hold on
        errorbar(mprA_0_30(:,1),mean(mprA_0_30(:,2:4),2)/mean(mprA0),std(mprA_0_30(:,2:4),'',2)/mean(mprA0),'k-^'); hold on
        errorbar(mprA_0_15(:,1),mean(mprA_0_15(:,2:4),2)/mean(mprA0),std(mprA_0_15(:,2:4),'',2)/mean(mprA0),'k--^')
        errorbar(mprA_30_0(:,1)+120,mean(mprA_30_0(:,2:4),2)/mean(mprA0),std(mprA_30_0(:,2:4),'',2)/mean(mprA0),'r-^')
        errorbar(mprA_30_15(:,1)+120,mean(mprA_30_15(:,2:4),2)/mean(mprA0),std(mprA_30_15(:,2:4),'',2)/mean(mprA0),'r--^')
        ylabel('normalized time-course'); xlabel('time(min)')
        title('mprA','FontSize',16);
        legend('0-->0.03','0-->0.015','0.03-->0','0.030-->0.015')
        xlim([-10 250])
    subplot(2,3,2)
        errorbar(sigE_wt(:,1),sigE_wt(:,2)/sigE_wt(1,2),sigE_wt_std/sigE_wt(1,2),'k-*')
        hold on
        errorbar(wt_up_E(:,1),mean(wt_up_E(:,2:4),2)/mean(wt_up_E(1,2:4),2),std((wt_up_E(:,2:4)),0,2)/mean(wt_up_E(1,2:4),2),'k-*');
        errorbar(wt_down_E(:,1)+120,mean(wt_down_E(:,2:4),2)/mean(wt_up_E(1,2:4),2),std(wt_down_E(:,2:4),0,2)/mean(wt_up_E(1,2:4),2),'r-*');
        errorbar(sigE_0_30(:,1),mean(sigE_0_30(:,2:4),2)/mean(sigE0),std(sigE_0_30(:,2:4),'',2)/mean(sigE0),'k-^'); hold on
        errorbar(sigE_30_0(:,1)+120,mean(sigE_30_0(:,2:4),2)/mean(sigE0),std(sigE_30_0(:,2:4),'',2)/mean(sigE0),'r-^')
        errorbar(sigE_0_15(:,1),mean(sigE_0_15(:,2:4),2)/mean(sigE0),std(sigE_0_15(:,2:4),'',2)/mean(sigE0),'k--^')
        errorbar(sigE_30_15(:,1)+120,mean(sigE_30_15(:,2:4),2)/mean(sigE0),std(sigE_30_15(:,2:4),'',2)/mean(sigE0),'r--^')
        xlabel('time(min)')
        title('sigE','FontSize',16)
        xlim([-10 250])
    subplot(2,3,3)
        errorbar(sigB_wt(:,1),mean(sigB_wt(:,2:4),2)/mean(sigB_wt(1,2:4),2),std(sigB_wt(:,2:4)/mean(sigB_wt(1,2:4),2),'',2),'k-*');
        hold on
        errorbar(sigB_wt2(:,1),mean(sigB_wt2(:,2:4),2)/mean(sigB_wt2(1,2:4),2),std(sigB_wt2(:,2:4)/mean(sigB_wt2(1,2:4),2),'',2),'k-*')
        errorbar(sigB_wt3(:,1),mean(sigB_wt3(:,2:4),2)/mean(sigB_wt3(1,2:4),2),std(sigB_wt3(:,2:4)/mean(sigB_wt3(1,2:4),2),'',2),'k-*')
        errorbar(sigB_0_30(:,1),mean(sigB_0_30(:,2:4),2)/mean(sigB0),std(sigB_0_30(:,2:4),'',2)/mean(sigB0),'k-^'); hold on
        errorbar(sigB_0_15(:,1),mean(sigB_0_15(:,2:4),2)/mean(sigB0),std(sigB_0_15(:,2:4),'',2)/mean(sigB0),'k--^')
        errorbar(sigB_30_0(:,1)+120,mean(sigB_30_0(:,2:4),2)/mean(sigB0),std(sigB_30_0(:,2:4),'',2)/mean(sigB0),'r-^')
        errorbar(sigB_30_15(:,1)+120,mean(sigB_30_15(:,2:4),2)/mean(sigB0),std(sigB_30_15(:,2:4),'',2)/mean(sigB0),'r--^')
        xlabel('time(min)')
        title('sigB','FontSize',16)
        xlim([-10 250])
    subplot(2,3,4)
        errorbar(mprA_sds_inc(:,1),mean(mprA_sds_inc(:,2:4),2)/mean(mprA_sds_inc(1,2:4),2),std(mprA_sds_inc(:,2:4),'',2)/mean(mprA_sds_inc(1,2:4),2),'k-*')
        hold on
        errorbar(mprA_sds_dec(:,1),mean(mprA_sds_dec(:,2:4),2)/mean(mprA_sds_inc(1,2:4),2),std(mprA_sds_dec(:,2:4),'',2)/mean(mprA_sds_inc(1,2:4),2),'r-*')
        errorbar(mprA_sds_inc2(:,1),mean(mprA_sds_inc2(:,2:4),2)/mean(mprA0),std(mprA_sds_inc2(:,2:4),'',2)/mean(mprA0),'k--^'); hold on
        errorbar(mprA_sds_dec2(:,1),mean(mprA_sds_dec2(:,2:4),2)/mean(mprA0),std(mprA_sds_dec2(:,2:4),'',2)/mean(mprA0),'r--*')
% plotting time series end points (240 min) for up and down sds
        errorbar(0.03,mean(mprA_0_30(end,2:4),2)/mean(mprA0),std(mprA_0_30(end,2:4),'',2)/mean(mprA0),'ko','MarkerFaceColor','y')
        errorbar(0.015,mean(mprA_0_15(end,2:4),2)/mean(mprA0),std(mprA_0_15(end,2:4),'',2)/mean(mprA0),'ko','MarkerFaceColor','y')
        errorbar(0,mean(mprA_30_0(end,2:4),2)/mean(mprA0),std(mprA_30_0(end,2:4),'',2)/mean(mprA0),'rs','MarkerFaceColor','g')
        errorbar(0.015,mean(mprA_30_15(end,2:4),2)/mean(mprA0),std(mprA_30_15(end,2:4),'',2)/mean(mprA0),'rs','MarkerFaceColor','g')
        ylabel('normalized concentration course'); xlabel('SDS [%w/v]')
        xlim([-0.004 0.035])
        legend('SDS up','SDS down','SDS up time course (2hr)','SDS down time course (2hr)')
    subplot(2,3,5)
        errorbar(sigE_sds_inc(:,1),mean(sigE_sds_inc(:,2:4),2)/mean(sigE_sds_inc(1,2:4),2),std(sigE_sds_inc(:,2:4),'',2)/mean(sigE_sds_inc(1,2:4),2),'k-*')
        hold on
        errorbar(sigE_sds_dec(:,1),mean(sigE_sds_dec(:,2:4),2)/mean(sigE_sds_inc(1,2:4),2),std(sigE_sds_dec(:,2:4),'',2)/mean(sigE_sds_inc(1,2:4),2),'r-*')
        errorbar(sigE_sds_inc2(:,1),mean(sigE_sds_inc2(:,2:4),2)/mean(sigE_sds_inc2(1,2:4),2),std(sigE_sds_inc2(:,2:4),'',2)/mean(sigE_sds_inc2(1,2:4),2),'k-s'); hold on
        errorbar(sigE_sds_dec2(:,1),mean(sigE_sds_dec2(:,2:4),2)/mean(sigE_sds_inc2(1,2:4),2),std(sigE_sds_dec2(:,2:4),'',2)/mean(sigE_sds_inc2(1,2:4),2),'r-s')
        errorbar(sigE_sds_inc3(:,1),mean(sigE_sds_inc3(:,2:4),2)/mean(sigE0),std(sigE_sds_inc3(:,2:4),'',2)/mean(sigE0),'k--^'); hold on
        errorbar(sigE_sds_dec3(:,1),mean(sigE_sds_dec3(:,2:4),2)/mean(sigE0),std(sigE_sds_dec3(:,2:4),'',2)/mean(sigE0),'r--^')
        errorbar(0.03,mean(sigE_0_30(end,2:4),2)/mean(sigE0),std(sigE_0_30(end,2:4),'',2)/mean(sigE0),'ko','MarkerFaceColor','y')
        errorbar(0.015,mean(sigE_0_15(end,2:4),2)/mean(sigE0),std(sigE_0_15(end,2:4),'',2)/mean(sigE0),'ko','MarkerFaceColor','y')
        errorbar(0.015,mean(sigE_30_15(end,2:4),2)/mean(sigE0),std(sigE_30_15(end,2:4),'',2)/mean(sigE0),'rs','MarkerFaceColor','g')
        errorbar(0,mean(sigE_30_0(end,2:4),2)/mean(sigE0),std(sigE_30_0(end,2:4),'',2)/mean(sigE0),'rs','MarkerFaceColor','g')
        xlabel('SDS [%w/v]')
        xlim([-0.004 0.035])
    subplot(2,3,6)
        errorbar(sigB_sds_inc(:,1),mean(sigB_sds_inc(:,2:4),2)/mean(sigB_sds_inc(1,2:4),2),std(sigB_sds_inc(:,2:4),'',2)/mean(sigB_sds_inc(1,2:4),2),'k-*')
        hold on
        errorbar(sigB_sds_dec(:,1),mean(sigB_sds_dec(:,2:4),2)/mean(sigB_sds_inc(1,2:4),2),std(sigB_sds_dec(:,2:4),'',2)/mean(sigB_sds_inc(1,2:4),2),'r-*')
        errorbar(sigB_sds_inc2(:,1),mean(sigB_sds_inc2(:,2:4),2)/mean(sigB0),std(sigB_sds_inc2(:,2:4),'',2)/mean(sigB0),'k--^'); hold on
        errorbar(sigB_sds_dec2(:,1),mean(sigB_sds_dec2(:,2:4),2)/mean(sigB0),std(sigB_sds_dec2(:,2:4),'',2)/mean(sigB0),'r--^')
        errorbar(0.03,mean(sigB_0_30(end,2:4),2)/mean(sigB0),std(sigB_0_30(end,2:4),'',2)/mean(sigB0),'ko','MarkerFaceColor','y')
        errorbar(0,mean(sigB_30_0(end,2:4),2)/mean(sigB0),std(sigB_30_0(end,2:4),'',2)/mean(sigB0),'rs','MarkerFaceColor','g')
        errorbar(0.015,mean(sigB_0_15(end,2:4),2)/mean(sigB0),std(sigB_0_15(end,2:4),'',2)/mean(sigB0),'ko','MarkerFaceColor','y')
        errorbar(0.015,mean(sigB_30_15(end,2:4),2)/mean(sigB0),std(sigB_30_15(end,2:4),'',2)/mean(sigB0),'rs','MarkerFaceColor','g')
        xlabel('SDS [%w/v]')
        xlim([-0.004 0.035])
        
        
        
ydata = [mean(sigE_sds_inc(1,2:4),2)/mean(sigE_sds_inc(1,2:4) ,2) mean(sigE_sds_dec(end,2:4),2)/mean(sigE_sds_inc(1,2:4) ,2)
         mean(sigE_sds_inc2(1,2:4),2)/mean(sigE_sds_inc2(1,2:4) ,2) mean(sigE_sds_dec2(end,2:4),2)/mean(sigE_sds_inc2(1,2:4) ,2)
         mean(sigE_sds_inc3(1,2:4),2)/mean(sigE_sds_inc3(1,2:4) ,2) mean(sigE_sds_dec3(end,2:4),2)/mean(sigE_sds_inc3(1,2:4) ,2)];
err =   [std(sigE_sds_inc(1,2:4),'',2)/mean(sigE_sds_inc(1,2:4) ,2) std(sigE_sds_dec(end,2:4),'',2)/mean(sigE_sds_inc(1,2:4) ,2)
         std(sigE_sds_inc2(1,2:4),'',2)/mean(sigE_sds_inc2(1,2:4),2) std(sigE_sds_dec2(end,2:4),'',2)/mean(sigE_sds_inc2(1,2:4) ,2)
         std(sigE_sds_inc3(1,2:4),'',2)/mean(sigE_sds_inc3(1,2:4) ,2) std(sigE_sds_dec3(end,2:4),'',2)/mean(sigE_sds_inc3(1,2:4) ,2)];
figure;
bar(1:3, ydata, 'grouped')
set(gca,'XTicklabel',{'I','II','III'})
ylabel('\textit{sigE} (fold-change)')
hold on
ngroups = size(ydata, 1);
nbars = size(ydata, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, ydata(:,i), err(:,i), '.');
    er.Color = [0 0 0];
end
legend('unstressed cells','pre-stressed cells')