close all
d=load('SlopeSineFits.dat');

x=1:12;
for n=1:size(d,1)
    beta=d(n,2)+d(n,3)*sin(2*pi*x*d(n,4)+d(n,5));
    [bmax,MonMax(n)]=max(beta);
end

load VosTransects

vos=[];
for n=1:numel(d(:,1))
    idx=find(round([VosT.BackMop]) == d(n,1));
    if isempty(idx)
        vos(n)=NaN;
    else
        vos(n)=VosT(idx(1)).Slope;
    end
end
figure('position',[126          77        1131         700]);
ax1=subplot(4,1,1);
idx=find(d(:,2) < 0.01);d(idx,2)=NaN;
plot(d(:,1),movmedian(d(:,2),11,'omitnan'),'b-','linewidth',2,'DisplayName','Surveys');
[a,b]=linfit(d(:,1),movmedian(d(:,2),11,'omitnan'));
hold on;
plot(d(:,1),d(:,1)*a+b,'b--','linewidth',2,'DisplayName','Survey Linear Fit');
plot(d(:,1),movmedian(vos,11,'omitnan'),'r-','linewidth',2,'DisplayName','CoastSat');
[a,b]=linfit(d(:,1),movmedian(vos,11,'omitnan')');
plot(d(:,1),d(:,1)*a+b,'r--','linewidth',2,'DisplayName','CoastSat Linear Fit');
% plot(d(:,1),d(:,6),'g-','linewidth',1,'DisplayName','Beta Lower');
% plot(d(:,1),d(:,10),'r-','linewidth',1,'DisplayName','Beta Upper');
set(gca,'xlim',[500 920],'ylim',[0.01 0.09],'fontsize',12);
%xlabel('Mop Number');
ylabel('Mean Slope');grid on;
title({'Sine Wave Fit:  3-month Running Mean MSL-MHW Beach Face Slopes',...
    'Approx. 1km (11 Mop) Running Median Applied to Slope, Amplitude, and Period Data'},'fontsize',16)
legend('location','northwest','NumColumns',4)

ax2=subplot(4,1,2);
plot(d(:,1),movmedian(d(:,3),11),'b-','linewidth',2,'DisplayName','Beta MSL-MHW');
hold on;
[a,b]=linfit(d(:,1),movmedian(d(:,3),11,'omitnan'));
hold on;
plot(d(:,1),d(:,1)*a+b,'b--','linewidth',2,'DisplayName','Survey Linear Fit');
% plot(d(:,1),d(:,7),'g-','linewidth',1,'DisplayName','Beta Lower');
% plot(d(:,1),d(:,11),'r-','linewidth',1,'DisplayName','Beta Upper');
set(gca,'xlim',[500 920],'ylim',[0 0.06],'fontsize',12);
%xlabel('Mop Number');
ylabel('Amplitude');grid on;

ax3=subplot(4,1,3);
plot(d(:,1),movmedian(1./d(:,4),11),'b-','linewidth',2,'DisplayName','Beta MSL-MHW');
hold on;
% plot(d(:,1),1./d(:,8),'g-','linewidth',1,'DisplayName','Beta Lower');
% plot(d(:,1),1./d(:,12),'r-','linewidth',1,'DisplayName','Beta Upper');
set(gca,'xlim',[500 920],'ylim',[0 24],'fontsize',12);
%xlabel('Mop Number');
ylabel('Period (Months)');grid on;

ax4=subplot(4,1,4);
plot(d(:,1),MonMax,'b.','markersize',10,'DisplayName','Beta MSL-MHW');hold on;
% plot(d(:,1),movmin(d(:,5),5),'b:','linewidth',2,'DisplayName','Beta MSL-MHW');hold on;
% plot(d(:,1),movmax(d(:,5),5),'b-','linewidth',2,'DisplayName','Beta MSL-MHW');
hold on;
% plot(d(:,1),d(:,9),'g-','linewidth',1,'DisplayName','Beta Lower');
% plot(d(:,1),d(:,13),'r-','linewidth',1,'DisplayName','Beta Upper');
set(gca,'xlim',[500 920],'ylim',[.5 12.5],'fontsize',12);
xlabel('Mop Number');ylabel('Max Slope Month');grid on;

dy=0.05;
pos=get(ax4,'position');set(ax4,'position',[pos(1) pos(2)-dy pos(3) pos(4)+dy/2]);
pos=get(ax3,'position');set(ax3,'position',[pos(1) pos(2)-dy/1. pos(3) pos(4)+dy/2]);
pos=get(ax2,'position');set(ax2,'position',[pos(1) pos(2)-dy/1. pos(3) pos(4)+dy/2]);
pos=get(ax1,'position');set(ax1,'position',[pos(1) pos(2)-dy/1. pos(3) pos(4)+dy/2]);

makepng('SlopeSineFitsVsMop.png')