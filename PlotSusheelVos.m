
% load v1 B A S data based on VOS
fid=fopen('CAmopSlopeParameters.txt','r');
 data = textscan(fid,'%d %s %5.3f %6.4f %d');
 mop=data{1}; % mops match to Vos transects
 b=(data{3}); % mean slopes from Vos

 % load Susheel B A S
load VosTransects-filt.mat
VosT(end,:)=[];

%% Figure 1  Vos vs Susheel at each mop with data

% make two empty mop vs slope vectors
KV=nan(1,12000);SA=KV;
% place mops with slope info into vectors
SA(round([VosT.BackMop{:}]))=[VosT.B]'; % Sushell mean slopes vs mop number
KV(round(mop))=b; % mean slopes from Vos

% compare mean slopes
figure('position',[68   354   560   420]);
plot(KV,SA,'+'); hold on;
set(gca,'xlim',[0 .25],'ylim',[0 .25]);
grid on;
plot([0 .25],[0 .25],'k-');
% make lin fit to any combined slope data below 0.15
idx=find( KV < 0.15 &  SA < 0.15);
[a b]=linfit(KV(idx)',SA(idx)');
plot(KV(idx),a.*KV(idx)+b,'r-');
title(['Linear Fit a=' num2str(a) '  ; b=' num2str(b)],'fontsize',14)
xlabel('Vos B');ylabel('Susheel B');

 %% Figure 2  Susheel ratio of Amplitude to Mean SLope

 pos=get(gcf,'position');figure('position',[pos(1)+100 pos(2)-100 pos(3:4)]);
 plot([VosT.B],[VosT.A],'*')
 idx=find(~isnan([VosT.B]) & ~isnan([VosT.A]) & [VosT.B] > 0);
 [a b]=linfit([VosT.B(idx)],[VosT.A(idx)]);
 hold on
 plot([VosT.B(idx)],a.*[VosT.B(idx)]+b,'k--')
 hold on; grid on;
 ylabel('Susheel Amplitude A');xlabel('Susheel Mean Slope B');
 text(0,.2,{'   Linear Fit ratio of',...
     [' Amplitude to Mean Slope = ' num2str(a)]},...
     'fontsize',14)
 title('Susheel Ampltudes are ~1/5 * Mean Slope vs Guess at ~1/2 in v1')

 %% Figure 3  Susheel phase vs Mop number

 pos=get(gcf,'position');figure('position',[pos(1)+100 pos(2)-100 pos(3:4)]);
 plot([VosT.BackMop{:}],[VosT.S],'*');grid on;
 ylabel('Phase S');xlabel('Mop Number');

 %% Figure 4  Histogram of phases
 pos=get(gcf,'position');figure('position',[pos(1)+100 pos(2)-100 pos(3:4)]);
 histogram([VosT.S],0:5:365);ylabel('Count');xlabel('Phase S');
 title('No Clear Pattern of Winter vs Summer max slopes?','fontsize',14)

function [a,b]= linfit(x,y)

a=0;
b=0;

if(length(x) > 0 & length(y) > 0)

X = [ones(length(x),1) x];
ab = X\y;
b=ab(1);
a=ab(2);

end

end