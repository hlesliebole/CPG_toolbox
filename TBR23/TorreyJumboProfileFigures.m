% 
clearvars
MopNumber=585;
matfile=['M' num2str(MopNumber,'%5.5i') 'SM.mat'];
load(matfile,'SM');

% identify jumbos
jumbo=find(contains({SM.File},'umbo'));

% reduce local SM struct array to just jumbo surveys
SM=SM(jumbo);

figure;
hold on;xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88)');
title(['Mop ' num2str(MopNumber) ': All Mop Jumbo Transect Profiles'],'fontsize',18);
for n=1:size(SM,2)
    plot(SM(n).X1D,SM(n).Z1Dtransect,'-');hold on
end
set(gca,'xdir','reverse','fontsize',12);box on;grid on;
%makepng('MopJumboTransectProfiles.png')
%% 

nmop=0;
for MopNumber=585:585%520:683
    fprintf('%i\n',MopNumber)
    nmop=nmop+1;

matfile=['M' num2str(MopNumber,'%5.5i') 'SM.mat'];
load(matfile,'SM');

% identify jumbos
jumbo=find(contains({SM.File},'umbo'));

% reduce local SM struct array to just jumbo surveys
SM=SM(jumbo);

%  Get the beach year for each jumbo by adding 92 days to the 
%  surveys datenum to push Oct 1 into the next year
BeachYr=year(datetime([SM.Datenum]+92,'convertfrom','datenum'));

% Get the months of the surveys
SurvMon=month(datetime([SM.Datenum],'convertfrom','datenum'));

% Get the quarters of the surveys
SurvQtr=1+floor(SurvMon/4);

% Get the winter(ONDJFM=1)/summer(AMJJAS=2) season of the surveys
SurvSsn=ones*size(SurvQtr);SurvSsn(SurvQtr == 2 | SurvQtr == 3)=2;

%fprintf('Stepping through Oct 1, (yr-1)-Sep 30,(yr) beach years...\n')

% step though beach years and get mean annual profiles
byr=unique(BeachYr);
AnnZ=NaN(numel(byr),numel(SM(1).X1D));
ny=0;

for y=byr
    ny=ny+1;
    
    % mean monthly profiles
    MonthlyZ=NaN(12,numel(SM(1).X1D));
    for m=1:12
        idx=find(BeachYr ==  y & SurvMon == m);
        if ~isempty(idx) 
            if numel([SM(idx).Z1Dtransect]) > 0
        zm=mean(vertcat(SM(idx).Z1Dtransect),1,'omitnan');
        MonthlyZ(m,:)=mean(vertcat(SM(idx).Z1Dtransect),1,'omitnan');
            end
        end
    end
    
    % mean quarterly profiles
    QtrlyZ=NaN(4,numel(SM(1).X1D));
    for q=1:4
        n=1+3*(q-1);
        QtrlyZ(q,:)=mean(MonthlyZ(n:n+2,:),1,'omitnan');
    end
    
    % mean 2 season profiles
        SsnZ=NaN(2,numel(SM(1).X1D));
        SsnZ(1,:)=mean(QtrlyZ([1 4],:),1,'omitnan'); % fall-winter
        SsnZ(2,:)=mean(QtrlyZ([2 3],:),1,'omitnan'); % spring-summer
    
    % annual
        AnnZ(ny,:)=mean(SsnZ,1); % will = NaN if one or both seasons have no data
end
        
% global mean

GlobalZ=mean(AnnZ,1,'omitnan');

% number years with enough data for annual means

Gyrs=numel(find(sum(~isnan(AnnZ)') > 0));

GMP(nmop).Mop=MopNumber;
GMP(nmop).NumDataYrs=Gyrs;
GMP(nmop).X1D=SM(1).X1D;
GMP(nmop).Z1D=GlobalZ;

end

figure('position',[440   377   623   420]);
hold on;xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88)');
title(['Mop ' num2str(MopNumber) ': Annual Mean Jumbo Transect Profiles'],'fontsize',18);
for n=1:size(AnnZ,1)
    if numel(find(~isnan(AnnZ(n,:)))) > 0
    plot(SM(n).X1D,AnnZ(n,:),'-','DisplayName',num2str(byr(n)));hold on
    end
end
set(gca,'xdir','reverse','fontsize',12);box on;grid on;
lg=legend('location','eastoutside');title(lg,'Beach Year');
%makepng('MopAnnualJumboTransectProfiles.png')
%% 

figure;
hold on;xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88)');
title(['Mop ' num2str(MopNumber) ': Global Mean Jumbo Transect Profile'],'fontsize',18);
for n=1:1
    plot(GMP(n).X1D,GMP(n).Z1D,'-','DisplayName',num2str(byr(n)));hold on
end
set(gca,'xdir','reverse','fontsize',12);box on;grid on;
%makepng('MopGlobalJumboTransectProfiles.png')

%%  Alongshore Homogeneity North Torrey Transects

nmop=0;
MopStart=580;
MopEnd=589;
for MopNumber=MopStart:MopEnd % Torrey TBR23 reach
    fprintf('%i\n',MopNumber)
    nmop=nmop+1;

matfile=['M' num2str(MopNumber,'%5.5i') 'SM.mat'];
load(matfile,'SM');

% identify jumbos
jumbo=find(contains({SM.File},'umbo'));

% reduce local SM struct array to just jumbo surveys
SM=SM(jumbo);

%  Get the beach year for each jumbo by adding 92 days to the 
%  surveys datenum to push Oct 1 into the next year
BeachYr=year(datetime([SM.Datenum]+92,'convertfrom','datenum'));

% Get the months of the surveys
SurvMon=month(datetime([SM.Datenum],'convertfrom','datenum'));

% Get the quarters of the surveys
SurvQtr=1+floor(SurvMon/4);

% Get the winter(ONDJFM=1)/summer(AMJJAS=2) season of the surveys
SurvSsn=ones*size(SurvQtr);SurvSsn(SurvQtr == 2 | SurvQtr == 3)=2;

%fprintf('Stepping through Oct 1, (yr-1)-Sep 30,(yr) beach years...\n')

% step though beach years and get mean annual profiles
byr=unique(BeachYr);
AnnZ=NaN(numel(byr),numel(SM(1).X1D));
ny=0;

for y=byr
    ny=ny+1;
    
    % mean monthly profiles
    MonthlyZ=NaN(12,numel(SM(1).X1D));
    for m=1:12
        idx=find(BeachYr ==  y & SurvMon == m);
        if ~isempty(idx) 
            if numel([SM(idx).Z1Dtransect]) > 0
        zm=mean(vertcat(SM(idx).Z1Dtransect),1,'omitnan');
        MonthlyZ(m,:)=mean(vertcat(SM(idx).Z1Dtransect),1,'omitnan');
            end
        end
    end
    
    % mean quarterly profiles
    QtrlyZ=NaN(4,numel(SM(1).X1D));
    for q=1:4
        n=1+3*(q-1);
        QtrlyZ(q,:)=mean(MonthlyZ(n:n+2,:),1,'omitnan');
    end
    
    % mean 2 season profiles
        SsnZ=NaN(2,numel(SM(1).X1D));
        SsnZ(1,:)=mean(QtrlyZ([1 4],:),1,'omitnan'); % fall-winter
        SsnZ(2,:)=mean(QtrlyZ([2 3],:),1,'omitnan'); % spring-summer
    
    % annual
        AnnZ(ny,:)=mean(SsnZ,1); % will = NaN if one or both seasons have no data
end
        
% global mean

GlobalZ=mean(AnnZ,1,'omitnan');

% number years with enough data for annual means

Gyrs=numel(find(sum(~isnan(AnnZ)') > 0));

GMP(nmop).Mop=MopNumber;
GMP(nmop).NumDataYrs=Gyrs;
GMP(nmop).X1D=SM(1).X1D;
GMP(nmop).Z1D=GlobalZ;

end

% save global mean profile struct array
matfile=['Mops' num2str(MopStart) 'to' num2str(MopEnd) 'GMP.mat'];
save(matfile,'GMP')
%% 

% alongshore global profile differences

% elev range to consider in comparing global profiles
minZ=-8;
maxZ=2;
GMP(1).Xlag=0;
GMP(1).RMSD=0;
GMP(1).Xlength=numel(~isnan(GMP(1).Z1D(GMP(1).Z1D > minZ & GMP(1).Z1D < maxZ)));

for n=1:size(GMP,2)-1
    %GMP(n).Mop;
    
    if n == 1
        gyrs(n)=GMP(1).NumDataYrs;
        mop(n)=GMP(1).Mop;
        xlag(n)=0;
        rmsd(n)=0;
        xlength(n)=numel(~isnan(GMP(1).Z1D(GMP(1).Z1D > minZ & GMP(1).Z1D < maxZ))); 
    end

    gyrs(n+1)=GMP(n+1).NumDataYrs;
    mop(n+1)=GMP(n+1).Mop;
    xlength(n+1)=NaN;
    xlag(n+1)=NaN;
    rmsd(n+1)=NaN;
    
    if numel(find(~isnan(GMP(n).Z1D))) > 0
    if numel(find(~isnan(GMP(n+1).Z1D))) > 0
      [rmsd(n+1),xlag(n+1),xlength(n+1)]=GetRmsdGMP(GMP,GMP(n).Mop,GMP(n+1).Mop,minZ,maxZ);  
      GMP(n+1).Xlag=xlag(n+1);
      GMP(n+1).RMSD=rmsd(n+1);
      GMP(n+1).Xlength=xlength(n+1);
    end
    end
end

%% rmsd plot for torrey pines
figure('position',[191   457   767   312]);hold on;
xlabel('MOP Number');ylabel('RMSD (m) with Adjacent Mop');
idx=find(mop >= MopStart & mop <= MopEnd); 
plot(mop(idx),rmsd(idx),'.-','markersize',10,'DisplayName','RMSD');hold on;
set(gca,'fontsize',14)
title('North Torrey Alongshore Adjacent Global Profile RMSD','fontsize',16);
grid on;box;
% fill([550 555 555 550 550],[0 0 1 1 0],'g','FaceAlpha',.1,'DisplayName','Blacks North');
% fill([580 584 584 580 580],[0 0 1 1 0],'r','FaceAlpha',.1,'DisplayName','Ruby2D');
% legend('location','northwest')
%makepng('TorreyPinesAlongshoreGlobalRMSD.png')



%% plot number of years of annual profiles
figure('position',[191   457   767   312]);hold on;
xlabel('MOP Number');ylabel('Years of Data');
idx=find(mop >= MopStart & mop <= MopEnd); 
plot(mop(idx),gyrs(idx),'.-','markersize',10,'DisplayName','Number of Yrs');hold on;
set(gca,'fontsize',14,'xlim',[MopStart MopEnd])
title('North Torrey Pines : Years of Data for Global Profile Estimates','fontsize',18);
grid on;box;
% fill([550 555 555 550 550],[0 0 25 25 0],'g','FaceAlpha',.1,'DisplayName','Blacks North');
% fill([580 584 584 580 580],[0 0 25 25 0],'r','FaceAlpha',.1,'DisplayName','Ruby2D');
% fill([599 636 636 599 599],[0 0 25 25 0],'c','FaceAlpha',.1,'DisplayName','Del Mar');
% plot([636 636],[0 25],'k--','DisplayName','San Dieguito Inlet')
%legend('location','southwest');
%makepng('NorthCountyAlongshoreGlobalDataYears.png')

%% plot xshore spatial lags of adjacent mops
figure('position',[191   457   767   312]);hold on;
xlabel('MOP Number');ylabel('Adjacent Xshore Axis Lag (m)');
idx=find(mop >= MopStart & mop <= MopEnd); 
plot(mop(idx),xlag(idx),'.-','markersize',10,'DisplayName','Adjacent Xshore Axis Lag (m)');hold on;
set(gca,'fontsize',14,'xlim',[MopStart MopEnd])
title('North Torrey Pines : min RMSD Adjacent Xshore Axis Spatial Lag (m)','fontsize',18);
grid on;box;
% fill([550 555 555 550 550],[-50 -50 50 50 -50],'g','FaceAlpha',.1,'DisplayName','Blacks North');
% fill([580 584 584 580 580],[-50 -50 50 50 -50],'r','FaceAlpha',.1,'DisplayName','Ruby2D');
% fill([599 636 636 599 599],[-50 -50 50 50 -50],'c','FaceAlpha',.1,'DisplayName','Del Mar');
% plot([636 636],[-50 50],'k--','DisplayName','San Dieguito Inlet')
% legend('location','south');
%set(gca,'ytick',[-50:10:50])
%makepng('NorthCountyProfileSpatialLags.png')

%% plot xshore lengths of -8 to 2m global  profiles
figure('position',[191   457   767   312]);hold on;
xlabel('MOP Number');ylabel('Profile Length (m)');
idx=find(mop >= MopStart & mop <= MopEnd); 
plot(mop(idx),xlength(idx),'.-','markersize',10,'DisplayName','Profile Length');hold on;
set(gca,'fontsize',14,'xlim',[MopStart MopEnd])
title('North Torrey Pines : Xshore Global Profile Length (m) [z= -8 to 2m]','fontsize',18);
grid on;box;
% fill([550 555 555 550 550],[100 100 600 600 100],'g','FaceAlpha',.1,'DisplayName','Blacks North');
% fill([580 584 584 580 580],[100 100 600 600 100],'r','FaceAlpha',.1,'DisplayName','Ruby2D');
% fill([599 636 636 599 599],[100 100 600 600 100],'c','FaceAlpha',.1,'DisplayName','Del Mar');
% plot([636 636],[100 600],'k--','DisplayName','San Dieguito Inlet')
% legend('location','south');
% set(gca,'ytick',100:50:600)
%makepng('NorthCountyProfileLengths.png')
%%
figure
hold on;xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88)');
title('Alongshore Aligned Global Mean Jumbo Transect Profiles','fontsize',18);

set(gca,'xdir','reverse','fontsize',12);box on;grid on;

xlag=0;
for m=MopStart:MopEnd
   idx=find([GMP.Mop] == m);
   xlag=xlag+GMP(idx).Xlag;
   plot(GMP(idx).X1D+xlag,GMP(idx).Z1D,'-','DisplayName',[num2str(m) '  Xlag ' num2str(xlag) 'm']);hold on;
end
legend('location','northwest')
%% 

function [MinRmsd,xlag,xlength]=GetRmsdGMP(GMP,mop1,mop2,minZ,maxZ)

imop1=find([GMP.Mop] == mop1);
imop2=find([GMP.Mop] == mop2);
MinRmsd=Inf;
for dx=-50:50
xp1=round(GMP(imop1).X1D);
xp2=round(GMP(imop2).X1D+dx);
gp1=GMP(imop1).Z1D;
gp2=GMP(imop2).Z1D;
gp1(gp1 < minZ)=NaN;
gp1(gp1 > maxZ)=NaN;
gp2(gp2 < minZ)=NaN;
gp2(gp2 > maxZ)=NaN;
idx1=find(~isnan(gp1));
idx2=find(~isnan(gp2));

overlap=max([xp1(idx1(1)) xp2(idx2(1))]):min([xp1(idx1(end)) xp2(idx2(end))]);

gp1=gp1(ismember(xp1,overlap));
gp2=gp2(ismember(xp2,overlap)); 

diff=gp2-gp1;
rmsd=sqrt(sum(diff.^2,'omitnan')/numel(~isnan(diff)));
%fprintf('%6.1f  %8.6f  %i\n',dx,rmsd,numel(~isnan(diff)))
MinRmsd=min([MinRmsd rmsd]);
if rmsd == MinRmsd;xlag = dx;xlength=numel(overlap);end
end

end

