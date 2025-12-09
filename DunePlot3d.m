close all

figure;
zc=1.566; % mhhw contour

load M00010SM.mat
xl=[SM(1).X1D(1) SM(1).X1D(end)];
for n=1:size(SM,2)
    
    x0=intersections(xl,[zc zc],SM(n).X1D,SM(n).Z1Dtransect);
    if ~isempty(x0)
    x0=max(x0);
    plot3(SM(n).X1D-x0,SM(n).Datenum*ones(size(SM(n).X1D))*0,SM(n).Z1Dtransect,'k-');hold on;  
    end
end
plot3(xl,SM(n).Datenum*ones(size(xl)),zc*ones(size(xl)),'k--')
set(gca,'xdir','reverse');set(gca,'xlim',[-100 100],'zlim',[.77 7]);grid on;

load M00674SM.mat
xl=[SM(1).X1D(1) SM(1).X1D(end)];
for n=1:size(SM,2)
    
    x0=intersections(xl,[zc zc],SM(n).X1D,SM(n).Z1Dtransect);
    if ~isempty(x0)
    x0=max(x0);
    plot3(SM(n).X1D-x0,SM(n).Datenum*ones(size(SM(n).X1D))*0,SM(n).Z1Dtransect,'r-');hold on;  
    end
end
n=size(SM,2)-1;
x0=intersections(xl,[zc zc],SM(n).X1D,SM(n).Z1Dtransect);
x0=max(x0);
plot3(SM(n).X1D-x0,SM(n).Datenum*ones(size(SM(n).X1D))*0,SM(n).Z1Dtransect,'m-','linewidth',2);hold on;  

load M00031SM.mat%M00010SM.mat%M00540SM.mat%M00066SM.mat%
xl=[SM(1).X1D(1) SM(1).X1D(end)];
for n=1:size(SM,2)
    
    x0=intersections(xl,[zc zc],SM(n).X1D,SM(n).Z1Dtransect);
    if ~isempty(x0)
    x0=max(x0);
    plot3(SM(n).X1D-x0,SM(n).Datenum*ones(size(SM(n).X1D))*0-1,SM(n).Z1Dtransect,'k-','linewidth',2);hold on;  
    end
end
plot3(xl,SM(n).Datenum*ones(size(xl)),zc*ones(size(xl)),'k--')
set(gca,'xdir','reverse');set(gca,'xlim',[-100 100],'zlim',[.77 7]);grid on;
% 
% figure;
% load M00674SM.mat
% 
% for n=1:size(SM,2)
%     plot(SM(n).X1D,SM(n).Z1Dtransect,'-','color',[.8 .8 .8]);hold on;
% end
% % plot(SM(end-1).X1D-12,SM(end-1).Z1Dmean,'k-','linewidth',2);hold on;  
% %set(gca,'xdir','reverse');set(gca,'xlim',[-100 100],'ylim',[.77 7]);grid on;
% 
% load M00674SM.mat
% %idx=find([SM.Datenum] > datenum(2019,8,28));
% idx=find([SM.Datenum] > datenum(2021,4,28));
% SM=SM(idx);
% 
% for n=1:size(SM,2)
%     plot(SM(n).X1D-12,SM(n).Z1Dtransect,'r-');hold on;
% end
% plot(SM(end-1).X1D-12,SM(end-1).Z1Dtransect,'k-','linewidth',2);hold on;  
% set(gca,'xdir','reverse');set(gca,'xlim',[-100 100],'ylim',[.77 7]);grid on;
% 
% load M00674SM.mat
% %yrs=year(datetime([SM.Datenum],'convertfrom','datenum'));
% %idx=find([SM.Datenum] < datenum(2018,10,9));
% idx=find([SM.Datenum] > datenum(2004,11,1) & [SM.Datenum] < datenum(2009,11,1));
% SM=SM(idx);
% %figure;
% for n=1:size(SM,2)
%     plot(SM(n).X1D,SM(n).Z1Dtransect,'b-');hold on;
% end
% set(gca,'xdir','reverse');set(gca,'xlim',[-100 100],'ylim',[.77 7]);grid on;
% 
% load M00674SM.mat
% %yrs=year(datetime([SM.Datenum],'convertfrom','datenum'));
% %idx=find([SM.Datenum] < datenum(2018,10,9));
% idx=find([SM.Datenum]  < datenum(2001,1,1));
% SM=SM(idx);
% %figure;
% for n=1:size(SM,2)
%     plot(SM(n).X1D,SM(n).Z1Dtransect,'k-','linewidth',2);hold on;
% end
% set(gca,'xdir','reverse');set(gca,'xlim',[-100 100],'ylim',[.77 7]);grid on;
% 
% %load M00010SM.mat
% load M00550SM.mat
% 
% for n=1:size(SM,2)
%     plot(SM(n).X1D-0,SM(n).Z1Dmean,'g-');hold on; 
%     %plot(SM(n).X1D+15,SM(n).Z1Dmean,'g-');hold on;   
% end
% set(gca,'xdir','reverse');set(gca,'xlim',[-100 100],'ylim',[.77 7]);grid on;