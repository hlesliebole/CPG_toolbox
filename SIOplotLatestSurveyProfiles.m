clearvars
close all
addpath /Users/William/Desktop/Mops
addpath /Users/William/Desktop/Mops/toolbox
%DefineMopPath

% loop throu SIO mops 507-514

for MopNumber=510:510%507:507
    

load(['M00' num2str(MopNumber,'%3.3i') 'SA.mat'])
ndx=find([SA.Datenum] ==  datenum(2022,11,2) |...
    [SA.Datenum] ==  datenum(2022,12,20) |...
    [SA.Datenum] ==  datenum(2023,1,3) |...
    [SA.Datenum] ==  datenum(2023,9,12) |...
    [SA.Datenum] > datenum(2023,10,11));
%ndx=find([SA.Datenum] > datenum(2017,9,1) & [SA.Datenum] > datenum(2022,12,1));
% fall surveys
% mon=month(datetime([SA.Datenum],'convertfrom','datenum'));
% ndx=find(mon > 6 & mon < 10);
stype=cell(numel(ndx),1);

all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};

close all
figure('position',[1          55        1412         742]);
m=0;
for n=ndx 
     m=m+1;
     [x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
     v(m)=sum(z1di(z1di > 2 & z1di < 3.5),'omitnan');
     t(m)=SA(n).Datenum;
     %fprintf('%s %i\n')
        if n == ndx(end)
           
            
         p(m)=plot(x1d,z1di,'r-',...
        'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
        'linewidth',4);hold on;
    
        elseif n == ndx(end-1)
           
            %[x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
         p(m)=plot(x1d,z1di,'r--',...
        'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
        'linewidth',4);hold on;

         elseif n == ndx(end-2)
           
            %[x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
         p(m)=plot(x1d,z1di,'r:',...
        'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
        'linewidth',4);hold on;

         elseif n == ndx(1)  
            %[x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
         p(m)=plot(x1d,z1di,'b:',...
        'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
        'linewidth',4);hold on;
         elseif n == ndx(2)  
            %[x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
         p(m)=plot(x1d,z1di,'b--',...
        'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
        'linewidth',4);hold on;
         elseif n == ndx(3)  
            %[x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
         p(m)=plot(x1d,z1di,'b-',...
        'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
        'linewidth',4);hold on;
        else
        %[x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
         p(m)=plot(x1d,z1di,'LineStyle','-',...
        'Marker',all_marks{mod(1+round(11*rand),13)},...
        'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
        'linewidth',2);hold on;
        end
end
  m=m+1;
  [x1d,z1di]=GetNonGriddedProfile(MopNumber,2);
  eln=plot(x1d,z1di,'k-','linewidth',4,'displayname','1998 El Nino');
%yl=get(gca,'ylim');set(gca,'ylim',[0 yl(2)]);
 
grid on;
set(gca,'fontsize',12);xlabel('Distance From Mop Back Beach Point (m)');ylabel('Elevation (m, NAVD88)')
set(gca,'ylim',[0 3.5],'xdir','reverse','fontsize',14);

% title([{['Mop ' num2str(MopNumber)]},{['Surveys Since ' datestr(SA(ndx(2)).Datenum)]}],...
%    'fontsize',16);
title([{['Mop ' num2str(MopNumber)]},{'2022-23 vs 2023-24 Profile Retreat' }],...
   'fontsize',16);
xl=get(gca,'xlim');%xl(2)=90;set(gca,'xlim',xl);
plot(xl,[2.119 2.119],'k--');text(xl(2),2.26,' HAT','fontsize',14);
plot(xl,[1.566 1.566],'k--');text(xl(2),1.7,' MHHW','fontsize',14);
plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',14);
plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',14);
plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',14);  
%legend(p,'location','eastoutside','fontsize',14)
legend([eln p],'location','eastoutside','fontsize',12,'numcolumns',1);%ceil(size(SA,2)/60))
makepng(['SioMop' num2str(MopNumber) 'Profiles.png'])
  
end
% 
%   % make a mean profile for mops 507-510
%   
%   
% % code to process the oceanside citizen surveys in MOP CPG files
% %
% %  Survey regions is approx Mops 866 to 926
% 
% 
% 
% MopNumber=866;
% 
% load(['M00' num2str(MopNumber,'%2.2i') 'SA.mat'])
% ndx=find([SA.Datenum] > datenum(2020,10,1));
% stype=cell(numel(ndx),1);
% 
% all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};
% 
% figure('position',[147         250        1092         520]);
% m=0;
% for n=ndx 
%      m=m+1;
%         if n == ndx(end)
%            
%             [x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
%          p(m)=plot(x1d,z1di,'m*-',...
%         'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
%         'linewidth',2);hold on;
%         else
%         [x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
%          p(m)=plot(x1d,z1di,'LineStyle','-',...
%         'Marker',all_marks{mod(1+round(11*rand),13)},...
%         'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
%         'linewidth',1);hold on;
%         end
% end
%     
% %yl=get(gca,'ylim');set(gca,'ylim',[0 yl(2)]);
%  
% grid on;
% set(gca,'fontsize',12);xlabel('MOP Number');ylabel('Elevation (m, NAVD88)')
% set(gca,'ylim',[-1 4],'xdir','reverse');
% 
% title([{['Mop ' num2str(MopNumber)]},{'Surveys Since Oct 2020'}],...
%    'fontsize',16);
% xl=get(gca,'xlim');xl(2)=100;set(gca,'xlim',xl);
% plot(xl,[2.119 2.119],'k--');text(xl(2),2.26,' HAT','fontsize',14);
% plot(xl,[1.566 1.566],'k--');text(xl(2),1.7,' MHHW','fontsize',14);
% plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',14);
% plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',14);
% plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',14);  
% legend(p,'location','eastoutside')
% makepng(['OceansideMop' num2str(MopNumber) 'Profiles.png'])
% 
% 
% 
% 
% 
% figure('position',[48         213        1246         518]);
%  %----- add cobble sightings
%  
% matfile=['M' num2str(MopNumber,'%5.5i') 'SA.mat'];
% load(matfile,'SA');
%  %  load Mop Transect Info
% load('MopTableUTM.mat','Mop');
% 
% % divide mop area into 20 mop subtransects at 1m xshore resolution,
% %  with an extra 100m of back beach for each
% [x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,20,[-100 0]);
% 
% idx=find(vertcat(SA.Class) > 1);
% Xutm=vertcat(SA.X);Yutm=vertcat(SA.Y);
% Xutm=Xutm(idx);Yutm=Yutm(idx);
% Z=vertcat(SA.Z);Z=Z(idx);
% dt=[];
% for n=1:size(SA,2)
%     dt=[dt' SA(n).Datenum*ones(size(SA(n).Z))']';
% end
% dt=dt(idx);
% 
% [dp,NearIdx]=...
%     pdist2([yst(:),xst(:)],[double(Yutm),double(Xutm)],'euclidean','smallest',1);
% 
% [row,col] = ind2sub(size(xst),NearIdx);
% 
% hold on;pc=plot(x1d(col),Z,'m.','DisplayName','Past ATV Cobble Sightings');
% 
% %--------------
% 
% load(['M' num2str(MopNumber,'%5.5i') 'SM.mat' ],'SM');
% 
% %figure;
% nn=0;
% ng=0;
% nskip=0;
% for n=1:nr+5
%     if n == nr+3
%         m=find([SM.Datenum] == datenum(2005,9,19));
%         SM(m).Z1Dtransect(SM(m).Z1Dtransect < 0)=NaN;
%     elseif n == nr+2
%         %m=find([SM.Datenum] == datenum(2021,9,22));
%         m=find([SM.Datenum] == datenum(2011,10,10));
%         SM(m).Z1Dtransect(SM(m).Z1Dtransect < 0)=NaN;
%     elseif n == nr+5
%         %m=find([SM.Datenum] == datenum(1997,10,12));
%         m=find([SM.Datenum] == datenum(2022,5,7));
%         SM(m).Z1Dtransect(SM(m).Z1Dtransect < 0)=NaN;
%     elseif n == nr+4
%         %m=find([SM.Datenum] == datenum(1997,10,12));
%         m=find([SM.Datenum] == datenum(1998,4,8));
%         SM(m).Z1Dtransect(SM(m).Z1Dtransect < 0)=NaN;
%     elseif n == nr+1
%         m=find([SM.Datenum] == datenum(2022,4,15));
%         SM(m).Z1Dtransect(SM(m).Z1Dtransect < 0)=NaN;
%     else
%        m=size(SM,2)+1-n-nskip;
%        while strcmpi(SM(m).Source,'iG8wheel') == 0
%            m=m-1;
%            nskip=nskip+1;
%        end
%     end
% %     if n==1;SM(m).Datenum=datenum(2021,10,11);end
% %     if n==2;SM(m).Datenum=datenum(2021,10,12);end
% %     if n==3;SM(m).Datenum=datenum(2021,10,13);end
%     
%     if ~isnan(min(SM(m).Z1Dmean)) 
%         
%     nn=nn+1;
%     z=SM(m).Z1Dtransect;z(z <-1.0)=NaN;
%     
%     %if numel(find(~isnan(z))) > 0
%         
%     xMSL=intersections([SM(m).X1D(1) SM(m).X1D(end)],[0.774 0.774],SM(m).X1D,SM(m).Z1Dtransect);
%     if isempty(xMSL);xMSL=999.9;end
%     if n == 1
%         p(nn)=plot(SM(m).X1D,z,'m-','linewidth',3,'DisplayName',...
%         [datestr(SM(m).Datenum,'mm/dd/yy') ' xMSL=' num2str(xMSL(end),'%4.1f') 'm ' SM(m).Source]);hold on;
%      
%     elseif n > nr 
%       if n == nr+2
%         p(nn)=plot(SM(m).X1D,z,'k:','linewidth',2,'DisplayName',...
%         [datestr(SM(m).Datenum,'mm/dd/yy') ' xMSL=' num2str(xMSL(end),'%4.1f') 'm ' SM(m).Source]);hold on;
%     elseif n == nr+5
%         p(nn)=plot(SM(m).X1D,z,'g-','linewidth',2,'DisplayName',...
%         [datestr(SM(m).Datenum,'mm/dd/yy') ' xMSL=' num2str(xMSL(end),'%4.1f') 'm South Swell Max xMSL']);hold on;
%       elseif n == nr+4
%         p(nn)=plot(SM(m).X1D,z,'r-','linewidth',2,'DisplayName',...
%         [datestr(SM(m).Datenum,'mm/dd/yy') ' xMSL=' num2str(xMSL(end),'%4.1f') 'm NASA/NOAA ATM LiDAR']);hold on;
%     elseif n == nr+3
%         p(nn)=plot(SM(m).X1D,z,'k-','linewidth',2,'DisplayName',...
%         [datestr(SM(m).Datenum,'mm/dd/yy') ' xMSL=' num2str(xMSL(end),'%4.1f') 'm ' SM(m).Source]);hold on;
%       else  
%         p(nn)=plot(SM(m).X1D,z,'b-','linewidth',2,'DisplayName',...
%         [datestr(SM(m).Datenum,'mm/dd/yy') ' xMSL=' num2str(xMSL(end),'%4.1f') 'm ' SM(m).Source]);hold on;  
%       end 
%     else
%         ng=ng+1;
%         p(nn)=plot(SM(m).X1D,z,'-','color',[.8 .8 .8]/ng,'linewidth',2,'DisplayName',...
%         [datestr(SM(m).Datenum,'mm/dd/yy') ' xMSL=' num2str(xMSL(end),'%4.1f') 'm ' SM(m).Source]);hold on;
%     end
%     end
%     %end
% end
% 
% % best historical fit index
% % bfidx=279;
% % nn=nn+1;
% % z=SM(bfidx).Z1Dmean;z(z < -0.7)=NaN;
% % p(nn)=plot(SM(bfidx).X1D,z,'k:','linewidth',2,'DisplayName',...
% %         [datestr(SM(bfidx).Datenum,'mm/dd/yy') ' Best Past Fit z=0.5-1.0m']);
%     
% xl=get(gca,'xlim');
% set(gca,'xlim',[0 xl(2)]);
% %yl=get(gca,'ylim');set(gca,'ylim',[0 yl(2)]);
% 
% plot(xl,[2.119 2.119],'k--');text(xl(2),2.26,' HAT','fontsize',14);
% plot(xl,[1.566 1.566],'k--');text(xl(2),1.7,' MHHW','fontsize',14);
% plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',14);
% plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',14);
% plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',14);
% %ps=plot(77,-0.31,'k.','markersize',20,'DisplayName','Paros');
% set(gca,'xdir','reverse','fontsize',14);grid on;box on;
% legend([p pc],'location','eastoutside');
% title(['MOP ' num2str(MopNumber) ' Transect Profiles']);
% xlabel('Cross-shore Distance (m)');
% ylabel('Elevation (m, NAVD88)');
% %set(gca,'ylim',[-1 4]);
% 
% if nr > 3
% makepng(['MOP' num2str(MopNumber) 'GBR22profiles.png'])
% else
% makepng(['MOP' num2str(MopNumber) 'GBR22profiles.png'])  
% end
% end  