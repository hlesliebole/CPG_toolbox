%ndx=find(strcmp({SG.Source},'RTKdrone'));
close all
clear all
MopNumber=582;
%load M00582SG.mat
load(['M' num2str(MopNumber,'%5.5i') 'SG.mat'],'SG');
mn=month(datetime([SA.Datenum],'convertfrom','datenum')); % month of surveys

% idx=find(strcmp({SG.Source},'AtvMR'))
% SG(idx)=[];
% idx=find(strcmp({SG.Source},'Gps'))
% SG=SG(idx);

% get global x,y limits for entire gridded data set
xmin=[];xmax=[];ymin=[];ymax=[];
for n=1:size(SG,2)
 np(n)=numel(SG(n).Z);
 xmin=min([xmin [SG(n).X]']);xmax=max([xmax [SG(n).X]']);
 ymin=min([ymin [SG(n).Y]']);ymax=max([ymax [SG(n).Y]']);
end

% x,y grid arrays
[X,Y]=meshgrid(xmin:xmax,ymin:ymax);
%zg=ng;
% zgmin=ng+1000;
% zgmax=ng-1000;

% number of surveys for each grid point
ng=X*0; % initialize count grid
for n=1:size(SG,2) % loop through surveys
idx=sub2ind(size(X),SG(n).Y-ymin+1,SG(n).X-xmin+1); %points with data
ng(idx)=ng(idx)+1; % add to counts
end

ndx=find(ng(:) > 0);% find grid points with at least 1 data points

zg=NaN(size(X));% initialize processed data grid
zgd=NaN([numel(ndx) size(SG,2)]);% create 2d array of 1D data-only points vs survey number
ngl=ng(ndx);% 2d array of counts for 1D data-only points vs survey number

for n=1:size(SG,2) % loop though surveys
    idx=sub2ind(size(X),SG(n).Y-ymin+1,SG(n).X-xmin+1);
    zg(idx)=SG(n).Z;
    zgd(:,n)=zg(ndx);
    zg=zg*NaN;
end
zg(ndx)=median(zgd','omitnan');
surf(zg);hold on;
zg(ndx)=min(zgd');
surf(zg);shading flat;
zg(ndx)=max(zgd');
surf(zg);shading flat;

% find the minimum sigma deviation for each gridded point
%  using a moving survey window size of "swin" consecutive 
%  surveys WITH DATA at a grid point.
swin=10; % survey window size
zgdsig=zgd*NaN; % initialize min sigma array as NaNs

for n=1:size(zgd,1) % loop through grid points with data
    zt=zgd(n,:);
    i=find(~isnan(zt));% surveys with data at this grid point
    if numel(i) >= swin % use desired running survey window if possible 
        win=swin;
    else
        win=numel(i); % otherwise use what there is
    end
    ztg=zt(i);ztsig=ztg*0+Inf; % reduce to just valid survey data time series
     for j=1:numel(i)-win+1 % loop through time series with window
         jw=j:j+win-1;
         % get sigma value for each data point
         sd=std(ztg(jw));mn=mean(ztg(jw));sig=abs((ztg(jw)-mn)./sd);
         % retain minimum sigma derived for each data point
         ztsig(jw)=min([ztsig(jw)' sig']');
     end   
     zgdsig(n,i)=ztsig; % save minimum sigmas for this grid point
end

% define median elevation surface for each grid point and
%  survey, using only the data points that have sigma deviations < 1,
%  and using a moving survey window.  In this case the survey window
%  is not based on the surveys with valid data, but the actual 
%  sequential surveys whether or not they actually have data at a
%  particular grid point.

swin=10; % survey moving window size
zgmed=zgd; % initialize median z array with original z data
zgmed(zgdsig(:) > 1)=NaN;% reduce to only grid-survey points with sigma < 1;
zgmed=movmedian(zgmed,swin,2,'omitnan');
%zg(ndx)=zgmed(:,10);surf(zg)

% now load 1m avg Survay SA struct array for this Mop 
load(['M' num2str(MopNumber,'%5.5i') 'SA.mat'],'SA');

% convert SA x,y values to grid indices

% for each survey in SA, calculate z deviations of SA points from the 
%  the corresponding median grid.
jdx=find(ismember({SA.File},{SG.File})); % index of SA files in SG
zmed=NaN(size(X));
zmed(ndx)=median(zgd','omitnan');
for n=jdx
    ism=find(strcmp({SG.File},SA(n).File));
    if ~isempty(ism)
%         zmed=NaN(size(X));
%         zmed(ndx)=zgmed(:,ism); % median grid for this survey
        idx=sub2ind(size(X),round(SA(n).Y)-ymin+1,round(SA(n).X)-xmin+1);
        dz=SA(n).Z-zmed(idx);SA(n).QC=dz;%/std(dz,'omitnan');
    else
        fprintf('No match: %s\n',SA(n).File)
    end
end

mm=24;
nn=jdx(mm);
zg=NaN(size(X));
figure;zg(ndx)=zgmed(:,mm);surf(X,Y,zg);hold on;plot3(vertcat(SA(nn).X),...
    vertcat(SA(nn).Y),vertcat(SA(nn).Z),'m.');
figure;plot(SA(nn).Z-SA(nn).QC,SA(nn).QC,'k.');

figure;zg(ndx)=zgmed(:,1);surf(X,Y,zg);hold on;plot3(vertcat(SA(jdx).X),...
    vertcat(SA(jdx).Y),vertcat(SA(jdx).QC),'m.');
figure;histogram(vertcat(SA(jdx).Z))
% save sigma deviations as new SA.QC=sigma struct array elements




%           
%  
% for n=1:size(zgd,2)
%     zt=zgd(:,n);
%     for l=1:1
%     %ibad=isoutlier(zt,'quartiles');zt(ibad)=NaN;
%     %plot(zdm(ibad,n),zdf(ibad,n),'g.');
%     ibad=isoutlier(zt,'mean');zt(ibad)=NaN;
%     %plot(zdm(ibad,n),zdf(ibad,n),'m.');  
%     ibad=isoutlier(zt,'median');zt(ibad)=NaN;
%     %plot(zdm(ibad,n),zdf(ibad,n),'r.');
%     end
%     zgd(:,n)=zt;
% end
%  
% %close all
% figure
% zgb=NaN([size(zgd,1) 1]); % initialize best value array
% zdiff=zgb; % initial diff from median array
% for n=1:size(zgd,1) % loop through 1d array of grid points with data
%     zmed=median(zgd(n,:),'omitnan'); % median value
%     sz=abs(zgd(n,:)-zmed); % diffs from median
%     %sz=sort(zgd(n,~isnan(zgd(n,:))));
%     %sz=sort(sz); % sorted diffs in ascending order
%     if numel(sz(~isnan(sz))) > 0
%        [dmin,imin]=min(sz);
%        zgb(n)=zgd(n,imin);
%        zdiff(n)=dmin;
%     end
% end
% zg(ndx)=zgb;%zdiff;%zgb;
% clf;surf(zg);shading flat;
% 
% % figure;plot(zgb,zdiff,'k.')
% % 
% % 
% % figure;plot(median(zgd','omitnan'),sum(~isnan(zgd')),'k.');
% % 
% % figure;plot(median(zgd','omitnan'),median(zgd','omitnan')-min(zgd'),'k.')
% % 
% % zdf=zgd;
% % zdm=zgd*NaN;
% % for n=1:size(zgd,2)
% %     zdf(:,n)=zgd(:,n)-zgb;
% %     zdm(:,n)=zgb;
% % end
% %  figure;plot(zdm(:),zdf(:),'k.')
% %  hold on;
% %  
% % for n=1:size(zgd,2)
% %     zt=zgd(:,n);
% %     for l=1:1
% %     ibad=isoutlier(zt,'quartiles');zt(ibad)=NaN;
% %     plot(zdm(ibad,n),zdf(ibad,n),'g.');
% %     ibad=isoutlier(zt,'mean');zt(ibad)=NaN;
% %     plot(zdm(ibad,n),zdf(ibad,n),'m.');  
% %     ibad=isoutlier(zt,'median');zt(ibad)=NaN;
% %     plot(zdm(ibad,n),zdf(ibad,n),'r.');
% %     end
% % end
% %  
% %  
% %  
% %  
% % figure;histogram(zgb-min(zgd')',[0:.01:10]);
% % figure;histogram(max(zgd')'-zgb,[0.5:.005:10]);
% % 
% % zbin=0:.005:10;
% % % --- minimum grid surface difference from median surface points
% % zdiff=zgb-min(zgd')';
% % % histogram of difference between min and median elev at each grid point
% % [np,~,binnum]=histcounts(zdiff,zbin);
% % [~,npk]=max(np); % find the peak bin 
% % % find the first histogram elev with 1 point or less past the peak
% % % this is the max allowable diff from median surface
% % zdmin=zbin(npk+min(find(np(npk+1:end) < 2)));
% % 
% % % --- maximum grid surface difference from median surface points
% % zdiff=max(zgd')'-zgb;
% % % histogram of difference between min and median elev at each grid point
% % [np,~,binnum]=histcounts(zdiff,zbin);
% % [~,npk]=max(np); % find the peak bin 
% % % find the first histogram elev with 1 point or less past the peak
% % % this is the max allowable diff from median surface
% % zdmax=zbin(npk+min(find(np(npk+1:end) < 2)));
% % 
% % 
% % 
