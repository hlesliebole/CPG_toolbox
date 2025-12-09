% code to compare last 2 jumbo survey profiles for a mop reach of
%  18 transects

clearvars
close all

addpath ../profiles

MopStart=630;
MopEnd=MopStart+34;% double 18

figure('position',[ 33          53        1381         739]);
ax1=axes('position',[.8 .05 .15 .9]);
xl=get(gca,'xlim');
hold on;
load MopTableUTM.mat
MopSB=Mop(MopStart:2:MopEnd,:);
for n=1:18
plot([MopSB.BackLon(n) MopSB.OffLon(n)],[MopSB.BackLat(n) MopSB.OffLat(n)],'m-','linewidth',2)
plot(MopSB.OffLon(n)-0.005,MopSB.OffLat(n),'k.','markersize',1)
text(MopSB.OffLon(n),MopSB.OffLat(n),MopSB.Name{n},...
    'horizontalalign','right','color','c','fontsize',15,'fontweight','bold')
end
plot_google_map('MapType','satellite')
axis off

MopName=MopSB.Name{18};MopNum1=str2num(MopName(3:5));

iii=0;
for nn=1:3
    for mm=1:6
        iii=iii+1;
        %MopName=MopSB.Name{19-iii};MopNum=str2num(MopName(3:5));
        MopName=MopSB.Name{iii};MopNum=str2num(MopName(3:5));
        axs(iii)=axes('position',[.06+(nn-1)*.23 .055+(mm-1)*.152 .21 .145],'box','on');
        set(axs(iii),'xlim',[-20 600],'xdir','reverse','ylim',[-2.2 2.2],...
            'ytick',-2:.5:2,'fontsize',12);
        grid on;

        yyaxis right;
        set(axs(iii),'xlim',[-20 600],'xdir','reverse','ylim',[-12 6.0],...
            'ytick',-10:2:4,'fontsize',12,'ycol','r');
        hold on;
        %MopNum=MopNum1+1-(nn-1)*6-(7-mm);%499+iii;
           
        if iii == 3% || iii == 9
            yyaxis left 

            ylabel('Elevation Change from Previous Jumbo (m)','fontsize',14);
            
        end

         if iii == 15% || iii == 9
            yyaxis right 

            ylabel('Elevation (m, NAVD88)','fontsize',14);
            
        end

     if iii > 6 & iii < 13
         yyaxis right
         set(axs(iii),'yticklabels',' ')
         yyaxis left
         set(axs(iii),'yticklabels',' ')

     end
    
    if iii == 7
        xlabel('Cross-Shore Distance (m)')
    else
         if (iii == 1 | iii == 7 | iii == 13)
         else          
            set(axs(iii),'xticklabels',' ')
         end
    end

        if iii > 20
            text(550,-14,['MOP ' num2str(MopNum)],'fontsize',14,'fontweight','bold')
        else
            yyaxis right
            text(500,2.6,['MOP ' num2str(MopNum)],'fontsize',16,'fontweight','bold')
            yyaxis left
        end

        if iii == 3
            
xl=get(gca,'xlim');%xl(2)=90;set(gca,'xlim',xl);
yyaxis left
plot(xl,[0 0],'k-');
yyaxis right
plot(xl,[.774 .774],'r-');text(xl(2),1.,' MSL','color','r','fontsize',14);
        else
            xl=get(gca,'xlim');
yyaxis left
plot(xl,[0 0],'k-');
yyaxis right
plot(xl,[.774 .774],'r-');
        end
            
load(['M00' num2str(MopNum,'%3.3i') 'SA.mat'])

hold on;
idx=find(contains({SA.File},'umbo') | contains({SA.File},'etski') );

n=idx(end-2);

% now pass to the profile function using a modest YdistTol that should
%  pull in most jumbo data that is near the transect
NumSubTrans=1; % 1 will just return the main transect info; 100 would make
                % statistical profiles with complete use of any 1m res LiDAR
XgapTol=5; % patch any cross-shore profile gaps < 5m wide
YdistTol=25; % can go as large as 50m up- downcoast for a single Mop area

% get nearest point profiles

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA(n),NumSubTrans,XgapTol,YdistTol);

%[x1d,z1di]=GetNonGriddedProfile(MopNum,n);
yyaxis right
pl3=plot(X1Dcpg,Z1Dtrans,'g-',...
        'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
        'linewidth',3);
x1=X1Dcpg;z1=Z1Dtrans;


n=idx(end-1);

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA(n),NumSubTrans,XgapTol,YdistTol);

%[x1d,z1di]=GetNonGriddedProfile(MopNum,n);
yyaxis right
pl4=plot(X1Dcpg,Z1Dtrans,'r-',...
        'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
        'linewidth',3);
set(gca,'ycol','r')
x2=X1Dcpg;z2=Z1Dtrans;
yyaxis left;
pl5=plot(x1,z2-z1,'k.-','linewidth',2,'DisplayName','Elev. Change');
set(gca,'ycol','k')

if iii == 12
    lg=legend([pl3 pl4 pl5],'location','northoutside','numcolumns',4,'fontsize',16);
    lg.Position=[0.2697 0.9692 0.2469 0.0277];
end

    end
end

makepng(['Mop' num2str(MopStart) 'to' num2str(MopEnd) datestr(SA(idx(end)).Datenum,'YYYYmmDD') '.png'])

