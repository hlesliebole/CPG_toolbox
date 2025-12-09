function GMP=GetGMPmopRange(MopStart,MopEnd)

nmop=0;
% MopStart=580;
% MopEnd=589;
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

%% 

% alongshore global profile differences

% elev range to consider in comparing global profiles
minZ=-8;
maxZ=2;
GMP(1).Xlag=0;
GMP(1).RMSD=0;
GMP(1).Xlength=numel(~isnan(GMP(1).Z1D(GMP(1).Z1D > minZ & GMP(1).Z1D < maxZ)));

for n=1:size(GMP,2)-1
    
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

end

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

