
% imperial beach nourishment
%   Sept 2012
%  344K m^3 design(380K observed).  0.53mm grains size (0.25mm native).  
%  +3.6-4m navd88 pad elevation. 1.5km long. 58m wide (mean).
% 

addpath ..
addpath ../..
addpath ../profiles/

clear pl
n=0;
figure
%load M00040SA.mat
nm=0;
zt=[];
for MopNumber=42:51%60
    nm=nm+1;
%for MopNumber=649:658%60
load(['M' num2str(MopNumber,'%5.5i') 'SA.mat'],'SA')


dt=datetime([SA.Datenum],'convertfrom','datenum');
idx=find(dt > datetime(2011,8,1) & dt < datetime(2017,12,1));
idx=find(dt > datetime(2001,8,1) & dt < datetime(2024,12,1));
% imperial beach nourishment 5 years
idx=find(dt > datetime(2012,10,1) & dt < datetime(2013,12,1));

% iyr=[2011 2012 2012 2013 2013 2014 2014 2015 2015 2016 2016 2017 2017];
% imo=[10 3 10 4 11 5 10 2 10 4 10 4 9];
% idy=[12 19 15 12 20 16 20 18 28 21 12 10 12];
% 
% idx=find(ismember([SA.Datenum],datenum(iyr,imo,idy,0,0,0)));

% solana
%idx=find(dt > datetime(2022,8,1) & dt < datetime(2024,12,1));

SA=SA(idx);
% now pass to the profile function using a modest YdistTol that should
%  pull in most jumbo data that is near the transect
NumSubTrans=5; % 1 will just return the main transect info; 100 would make
                % statistical profiles with complete use of any 1m res LiDAR
XgapTol=15; % patch any cross-shore profile gaps < 5m wide
YdistTol=25; % can go as large as 50m up- downcoast for a single Mop area

% get nearest point profiles

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);
% (Note: It can take a few minutes to run depending on the number of 
%    NumSubTrans subtransects being considered.)

ix0=find(X1Dmop == 0);
for nn=1:size(Z1Dmedian,1)
    zt(nn,1:150,nm)=Z1Dmedian(nn,ix0:ix0+149);
end

if MopNumber == 51
    figure('position',[138         192        1093         527]);
    clear plt
    col=jet(numel(1:2:size(zt,1)));
    col(1,:)=[0 0 0];
    zt=mean(zt,3,'omitnan');
    nc=0;
    for nn=1:2:size(zt,1)
        nc=nc+1;
        if nn == 1
            p1=plot(0:149,zt(nn,:),'-','color',col(nc,:),'linewidth',3,'DisplayName',['Imperial Beach 2012 : ' char(string(Zdatetime(nn)))]);
            hold on;
        else
          %plt(nc)=plot(0:149,zt(nn,:),'.-','color',col(nc,:),'linewidth',3,'DisplayName',string(Zdatetime(nn)));
          hold on;  
        end
    end
    plot([0 150],[.774 .774],'k--','linewidth',2);
    text(15,.5,'MSL','fontsize',18)
    text(120,-.95,'Base of Pad','fontsize',18)
    text(45,3.25,'Terrace','fontsize',18)
    
    set(gca,'xdir','reverse','fontsize',16,'color',[.8 .8 .8],'xtick',[0:25:150],'linewidth',2);grid on;
    set(gca,'ylim',[-2 7])
    ylabel('Elevation (NAVD88, m)');
    yyaxis right
    set(gca,'xdir','reverse','fontsize',16,'color',[.8 .8 .8],'xtick',[0:25:150],'linewidth',2);grid on;
    set(gca,'ylim',[-2 7],'ycolor','k')
    xlabel('Cross-shore Distance from the back of the beach (m)');ylabel('Elevation (NAVD88, m)');
    title({'Past and Present Nourishment Fall Profiles Prior to Their First Winter'},'fontsize',22)
end

% figure;
% pcolor(X1Dcpg,Zdatetime,Z1Dmedian);shading flat;
% set(gca,'ydir','normal');datetick('y');
% BeachColorbar; % add beach friendly colormap/colorbar from the toolbox


%figure

% subplot(2,1,1)
% Znavd88=1.344;
% [BWmin,BWmax]=GetCpgProfileBeachWidths(Znavd88,X1Dcpg,Z1Dmedian);
% n=n+1;
% pl(n)=plot(Zdatetime,BWmax,'*-','displayname',num2str(MopNumber));grid on;hold on;
% 
% subplot(2,1,2)
% Z1=1.6;
% [BWmin,BWmax1]=GetCpgProfileBeachWidths(Z1,X1Dcpg,Z1Dmedian);
% Z2=1.;
% [BWmin0,BWmax2]=GetCpgProfileBeachWidths(Z2,X1Dcpg,Z1Dmedian);
% beta=(Z1-Z2)./(BWmax2-BWmax1);
% plot(Zdatetime,beta,'*-','displayname',num2str(MopNumber));grid on;hold on;

end
%legend(pl,'location','eastoutside')

%% ------------------------------------------------------

zt=[];
nm=0;
for MopNumber=652:655%649:658%60
    nm=nm+1;
load(['M' num2str(MopNumber,'%5.5i') 'SA.mat'],'SA')


dt=datetime([SA.Datenum],'convertfrom','datenum');
idx=find(dt > datetime(2011,8,1) & dt < datetime(2017,12,1));
idx=find(dt > datetime(2001,8,1) & dt < datetime(2024,12,1));
% imperial beach nourishment 5 years
idx=find(dt > datetime(2012,10,1) & dt < datetime(2017,11,1));

iyr=[2024];
imo=[10];
idy=[9];

idx=find(ismember([SA.Datenum],datenum(iyr,imo,idy,0,0,0)));

% solana
%idx=find(dt > datetime(2022,8,1) & dt < datetime(2024,12,1));

SA=SA(idx);
% now pass to the profile function using a modest YdistTol that should
%  pull in most jumbo data that is near the transect
NumSubTrans=5; % 1 will just return the main transect info; 100 would make
                % statistical profiles with complete use of any 1m res LiDAR
XgapTol=15; % patch any cross-shore profile gaps < 5m wide
YdistTol=25; % can go as large as 50m up- downcoast for a single Mop area

% get nearest point profiles

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);
% (Note: It can take a few minutes to run depending on the number of 
%    NumSubTrans subtransects being considered.)

ix0=find(X1Dcpg == 0);
for nn=1:size(Z1Dmedian,1)
    zt(nn,1:150,nm)=Z1Dmedian(nn,ix0:ix0+149);
end

if MopNumber == 655%658
    %figure('position',[138         192        1093         527]);
    %clear plt
    %col=jet(numel(1:2:size(zt,1)));
    col(nc+1,:)=[0 .8 0];%[1 0 1];
    zt=mean(zt,3,'omitnan');
    %nc=0;
    for nn=1:2:size(zt,1)
        nc=nc+1;
        if nn == 1
            p3=plot([0:100 125],[zt(nn,1:101) .2],'-','color',col(nc,:),'linewidth',3,'DisplayName',['Solana Beach 2024  : ' char(string(Zdatetime(nn)))]);
            hold on;
        else
          plt(nc)=plot(0:149,zt(nn,:),'.-','color',col(nc,:),'linewidth',3,'DisplayName',string(Zdatetime(nn)));
          hold on;  
        end
    end
    % plot([0 150],[.774 .774],'k--','linewidth',2);
    % text(140,.95,'MSL','fontsize',18)
    % legend(plt,'location','northwest')
    % set(gca,'xdir','reverse','fontsize',16,'color',[.8 .8 .8],'xtick',[0:25:150],'linewidth',2);grid on;
    % xlabel('Cross-shore Distance (m)');ylabel('Elevation (NAVD88, m)');
    % title('Imperial Beach Nourishment Profile Evolution (1km Alongshore Average)')
end
% figure;
% pcolor(X1Dcpg,Zdatetime,Z1Dmedian);shading flat;
% set(gca,'ydir','normal');datetick('y');
% BeachColorbar; % add beach friendly colormap/colorbar from the toolbox


%figure

% subplot(2,1,1)
% Znavd88=1.344;
% [BWmin,BWmax]=GetCpgProfileBeachWidths(Znavd88,X1Dcpg,Z1Dmedian);
% n=n+1;
% pl(n)=plot(Zdatetime,BWmax,'*-','displayname',num2str(MopNumber));grid on;hold on;
% 
% subplot(2,1,2)
% Z1=1.6;
% [BWmin,BWmax1]=GetCpgProfileBeachWidths(Z1,X1Dcpg,Z1Dmedian);
% Z2=1.;
% [BWmin0,BWmax2]=GetCpgProfileBeachWidths(Z2,X1Dcpg,Z1Dmedian);
% beta=(Z1-Z2)./(BWmax2-BWmax1);
% plot(Zdatetime,beta,'*-','displayname',num2str(MopNumber));grid on;hold on;

end

%% ------------------------------------------------------

zt=[];
nm=0;
for MopNumber=722:723%649:658%60
    nm=nm+1;
load(['M' num2str(MopNumber,'%5.5i') 'SA.mat'],'SA')


dt=datetime([SA.Datenum],'convertfrom','datenum');
idx=find(dt > datetime(2011,8,1) & dt < datetime(2017,12,1));
idx=find(dt > datetime(2001,8,1) & dt < datetime(2024,12,1));
% imperial beach nourishment 5 years
idx=find(dt > datetime(2012,10,1) & dt < datetime(2017,11,1));

iyr=[2024];
imo=[10];
idy=[16];

idx=find(ismember([SA.Datenum],datenum(iyr,imo,idy,0,0,0)));

% solana
%idx=find(dt > datetime(2022,8,1) & dt < datetime(2024,12,1));

SA=SA(idx);
% now pass to the profile function using a modest YdistTol that should
%  pull in most jumbo data that is near the transect
NumSubTrans=5; % 1 will just return the main transect info; 100 would make
                % statistical profiles with complete use of any 1m res LiDAR
XgapTol=15; % patch any cross-shore profile gaps < 5m wide
YdistTol=25; % can go as large as 50m up- downcoast for a single Mop area

% get nearest point profiles

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);
% (Note: It can take a few minutes to run depending on the number of 
%    NumSubTrans subtransects being considered.)

ix0=find(X1Dcpg == 0);
for nn=1:size(Z1Dmedian,1)
    zt(nn,1:150,nm)=Z1Dmedian(nn,ix0:ix0+149);
end

if MopNumber == 722%658
    %figure('position',[138         192        1093         527]);
    %clear plt
    %col=jet(numel(1:2:size(zt,1)));
    col(nc+1,:)=[1 0 1];
    zt=mean(zt,3,'omitnan');
    %nc=0;
    for nn=1:2:size(zt,1)
        nc=nc+1;
        if nn == 1
            %plot([105 125],[.7763 .5],'.-','color',col(nc,:),'linewidth',3)
            p2=plot([0:105 125],[zt(nn,1:106) .3],'-','color',col(nc,:),'linewidth',3,'DisplayName',['Moonlight SB 2024  : ' char(string(Zdatetime(nn)))]);
            %plot(0:149,zt(nn,:),'-','color',col(nc,:),'linewidth',3,'DisplayName',['Moonlight SB 2024  : ' char(string(Zdatetime(nn)))]);
            hold on;
        else
          plt(nc)=plot(0:149,zt(nn,:),'.-','color',col(nc,:),'linewidth',3,'DisplayName',string(Zdatetime(nn)));
          plt(nc)=plot([105 125],[.7763 .65],'.-','color',col(nc,:),'linewidth',3,'DisplayName',string(Zdatetime(nn)));
          hold on;  
          hold on;  
        end
    end
    % plot([0 150],[.774 .774],'k--','linewidth',2);
    % text(140,.95,'MSL','fontsize',18)
    % legend(plt,'location','northwest')
    % set(gca,'xdir','reverse','fontsize',16,'color',[.8 .8 .8],'xtick',[0:25:150],'linewidth',2);grid on;
    % xlabel('Cross-shore Distance (m)');ylabel('Elevation (NAVD88, m)');
    % title('Imperial Beach Nourishment Profile Evolution (1km Alongshore Average)')
end
% figure;
% pcolor(X1Dcpg,Zdatetime,Z1Dmedian);shading flat;
% set(gca,'ydir','normal');datetick('y');
% BeachColorbar; % add beach friendly colormap/colorbar from the toolbox


%figure

% subplot(2,1,1)
% Znavd88=1.344;
% [BWmin,BWmax]=GetCpgProfileBeachWidths(Znavd88,X1Dcpg,Z1Dmedian);
% n=n+1;
% pl(n)=plot(Zdatetime,BWmax,'*-','displayname',num2str(MopNumber));grid on;hold on;
% 
% subplot(2,1,2)
% Z1=1.6;
% [BWmin,BWmax1]=GetCpgProfileBeachWidths(Z1,X1Dcpg,Z1Dmedian);
% Z2=1.;
% [BWmin0,BWmax2]=GetCpgProfileBeachWidths(Z2,X1Dcpg,Z1Dmedian);
% beta=(Z1-Z2)./(BWmax2-BWmax1);
% plot(Zdatetime,beta,'*-','displayname',num2str(MopNumber));grid on;hold on;

end

legend([p2 p3 p1],'location','northwest','fontsize',18)

set(gcf,'color','w')
set(gcf,'inverthardcopy','off')
makepng('IBvSolanavMoonlight.png')

