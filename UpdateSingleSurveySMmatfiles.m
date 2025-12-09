% Adds a Single New Survey to the corresponding Mop SM matfiles
%

fprintf('UpdateSingleSurveySMmatfiles.m \n')

% SingleDate (a survey datenum) and SingleSource (survey source name)
% must be defined prio to runnnig the script.  Designed to be
% called by UpdateSingleSurveySASGSM.m

load SurveyMasterList.mat

idx=find([Survey.Datenum] == SingleDate & contains({Survey.Source},SingleSource) == 1);

load('MopTableUTM.mat','Mop');  % Load "Mop" table array

if numel(idx) == 1
    
Mop1=Survey(idx).MopStart;
Mop2=Survey(idx).MopEnd;
%Mop1=585;Mop2=585;

for MopID=Mop1:Mop2

    fprintf('Mop %i \n',MopID)
    
% Get numeric mop number
if isnumeric(MopID);MopNumber=MopID;else;...
        MopNumber=find(contains(Mop.Name,MopID));end

% load selected mop's MopSurvey struct array

matfile=['M' num2str(MopNumber,'%5.5i') 'SG.mat'];
if exist(matfile,'file') == 2
load(matfile,'SG');

SGidx=find([SG.Datenum] == datenum(2021,10,7) & contains({SG.Source},'Gps') == 1);
if numel(SGidx) == 1
    
% load exisiting grid SM struct array as old array OM
matfile=['M' num2str(MopNumber,'%5.5i') 'SM.mat'];
if exist(matfile,'file') == 2
 load(matfile,'SM');
 OM=SM;
else
 clear OM
 OM(1).Datenum=0;
 OM(2).Source=' ';
end

% SM array should be one survey smaller than SG. If
%  not, stop.

if size(SG,2) - size(SM,2) == 1

% Initializing new SM struct array up to the
%  survey idx-1 to match SA array
SM=OM(1:SGidx-1);

% New survey date
SurvNum=SGidx;
  
      SurveyDatenum=SG(SurvNum).Datenum;
      fprintf('Adding Survey %i %s\n',...
          SurvNum,datestr(SurveyDatenum))
      
 %-----------------------------------------------------------------   
 % divide mop area into 20 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,20,[-100 0]);


%-----------------------------------------------------------------
%  Individual survey morpho parameters
%---------------------------------------------------------------


%-----------------------------------------------------------------
% ---- 1D parameter fields ----
%-----------------------------------------------------------------
% .X1D : xshore distance (m) from Mop back beach line
% .Z1Dtransect : gridded z interpolated on Mop transect from gridded data  
% .Z1Dmean : mean gridded z at xshore distance X1D
% .Z1Dmedian : mean gridded z at xshore distance X1D
% .Z1Dstd : standard deviation gridded z at xshore distance X1D
% .Z1Dmin : minimum gridded z at xshore distance X1D
% .Z1Dmax : minimum gridded z at xshore distance X1D

% loop through new SG gridded data

for sn=SurvNum
    
    SM(sn).Mopnum=SG(sn).Mopnum;
    SM(sn).Datenum=SG(sn).Datenum;
    SM(sn).Source=SG(sn).Source;
    SM(sn).File=SG(sn).File;
    SM(sn).UTMzone=SG(sn).UTMzone;
    SM(sn).X1D=x1d;
    
% reconstruct x,y grid points within gridded x,y area that 
%  have "no data" NaN's

if ~isempty(SG(sn).X)
xg=SG(sn).X;yg=SG(sn).Y;zg=SG(sn).Z;cg=SG(sn).Class;
%cmin=accumarray([xidx(:), yidx(:)], cg.',[],@nanmax);
nx=max(xg)-min(xg)+3;ny=max(yg)-min(yg)+3; % temp grid dimensions
tg=nan(ny,nx); % temp grid of Nans
idx=sub2ind([ny nx],yg-min(yg)+2,xg-min(xg)+2); % data point 1d indices
tg(idx)=zg; % add data to temp grid
[X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1);

% now 2d interpolate z values for the Mop transect points
zt=xt*NaN; % initialize transect elevation points 
zt(:) = interp2(X,Y,tg,xt,yt);
  SM(sn).Z1Dtransect=zt;

% now 2d interpolate z values for all the subtransect points
zst=xst*NaN; % initialize transect elevation points 
zst(:) = interp2(X,Y,tg,xst(:),yst(:));
% get mean, median, std, min-max z(xt) transect values for 51 subtransects.

  SM(sn).Z1Dmean=nanmean(zst,1);
  SM(sn).Z1Dmedian=nanmedian(zst,1);
  SM(sn).Z1Dstd=nanstd(zst,1);
  SM(sn).Z1Dmin=nanmin(zst,[],1);
  SM(sn).Z1Dmax=nanmax(zst,[],1);
  
  tg(idx)=cg; % add class data to temp grid
  % now 2d interpolate class values for all the subtransect points
  cst=xst*NaN; % initialize transect class points 
  cst(:) = interp2(X,Y,tg,xst(:),yst(:));
  % find largest class id (hardest substrate)
  SM(sn).Z1Dclass=nanmax(cst,[],1);
end
   
end
 
     
% add any existing surveys after the added one
if size(OM,2)+1 > size(SM,2)
    SM(SGidx+1:size(OM,2)+1)=OM(SGidx:end);
end

% Save MopGrid struct array in its own mat file
matfile=['M' num2str(MopNumber,'%5.5i') 'SM.mat'];
save(matfile,'SM');

else
    fprintf('SG %i is not 1 survey larger than SM %i. Skipping.\n',...
        size(SG,2),size(SM,2))
end

else
   fprintf('SG does not contain the survey. Skipping.\n') 
end

else
   fprintf('SG mat files does not exist. Skipping.\n') 
end 

end  % end mop loop

else
    fprintf('No Matching Survey for that date and source.\n')
end