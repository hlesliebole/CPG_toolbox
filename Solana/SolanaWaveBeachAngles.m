% code to plot mop sxy map for a mop reach

clearvars
%close all

load SolanaShoreboxMap.mat

figure('position',[ 33          53        1381         739]);
image(xgimg,ygimg,sbgimg);set(gca,'ydir','normal');hold on;
set(gca,'dataaspectratio',[1 1 1],'xdir','normal','fontsize',16);
set(gca,'ylim',[min(ygimg) max(ygimg)],'xlim',[min(xgimg) max(xgimg)])
pos=get(gca,'position');
%set(gca,'position',[.5 pos(2) pos(3) pos(4)],'fontsize',14,'linewidth',2)
xlabel('Alongshore Distance (m)','fontsize',18);
%ylabel('Cross-shore Distance (m)','fontsize',18);
%title('Solana Beach Nourishment Pad Evolution :  10 May vs 16 Apr 2024','fontsize',22)
xl=get(gca,'xlim');
hold on;
for n=3:31
%plot([MopSB.BackLon(n) MopSB.OffLon(n)],[MopSB.BackLat(n) MopSB.OffLat(n)],'m-','linewidth',2)
%text(MopSB.OffLat(n)+300,MopSB.OffLon(n),[ '   ' MopSB.Name{n}],'color','y','fontsize',16,'fontweight','bold')
end
for n=4:2:30
%plot([MopSB.BackLat(n) MopSB.OffLat(n)],[MopSB.BackLon(n) MopSB.OffLon(n)],'m-','linewidth',2)
% plot(MopSB.OffLon(n),MopSB.OffLat(n),'y.','markersize',15)
% text(MopSB.OffLon(n),MopSB.OffLat(n),[ '  ' MopSB.Name{n}],'color','y','fontsize',18,'fontweight','bold','rotation',90)
end

tf=text(1850,-100,{'Fletcher','  Cove'},'color','w','fontsize',22,'fontweight','bold');
set(gca,'ylim',[-200 300],'ytick',[])
title('')
pos=get(gca,'position');
set(gca,'position',[pos(1) pos(2)-.3 pos(3) pos(4)])



ax2=axes('position',[pos(1) pos(2)+.3 pos(3) pos(4)-.3]);
% Sxy axes
%ax2=axes('position',[pos(1) .42 pos(3) .53 ],'fontsize',14,'linewidth',2);
%set(ax2,'xlim',xl);hold on;box on


% for m=635:665
%     stn=['D' num2str(m,'%4.4i') ];
%      [wavetime,Dm]=GetMopDm(stn);
%      SxyAngle(m-634)=mean(rmoutliers(Dm(year(wavetime) < 2024)));
%      SxyAngleSdev(m-634)=std(rmoutliers(Dm(year(wavetime) < 2024)));
%      SxyAngleFeb(m-634)=mean(Dm(month(wavetime) == 2));
%      SxyAngleAug(m-634)=mean(Dm(month(wavetime) == 8));
% end
% 
% save SolanaSxyAngles.mat SxyAngleSdev SxyAngle SxyAngleAug SxyAngleFeb

load SolanaSxyAngles.mat

% plot(635:665,SxyAngle,'k*-','linewidth',2,'DisplayName','Mop Mean Sxy Angle');
% plot(635:665,SxyAngleFeb,'b*-','linewidth',2,'DisplayName','Mop Mean Sxy Angle Feb');
% plot(635:665,SxyAngleAug,'g*-','linewidth',2,'DisplayName','Mop Mean Sxy Angle Aug');
errorbar(635:665,SxyAngle,SxyAngleSdev,'o-','linewidth',2,'DisplayName','Mop 2000-2023 Bulk Sxy Angle @ 10m Contour');
hold on



load SolanaLidarShoreNormals.mat 

errorbar(638:664,MopLnorm(4:end-1),MopLnormSdev(4:end-1),'o-','linewidth',2,...
    'DisplayName','LiDAR 2017-2024 MHW Beach Face Normal')


load MopTableUTM.mat

plot(635:665,Mop.Normal(635:665),'o-','linewidth',2,'DisplayName','Mop Fixed Shore Normal');grid on;
hold on;

xlabel('Mop Number');ylabel('Direction');set(gca,'xlim',[632 667],'xtick',635:665,'fontsize',16)
title('LiDAR MHW Beach Face Orientation and MOP incident Bulk Sxy Angle Means & Standard Deviations ')
box on;
legend('location','south','fontsize',16)
makepng('SolanaWaveBeachAngles.png')
