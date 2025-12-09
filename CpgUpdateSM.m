function []=CpgUpdateSM(UpdatedMops)

% For the mop numbers in the vector UpdatedMops,
% regenerate SM  Morpho struct files from SG grid files

CpgDefineMopPath
%  load Mop Transect Info
load('MopTableUTM.mat','Mop');


 for MopNumber=UpdatedMops
  if MopNumber > 1
  SM=[];
    %fprintf('%i of %i\n',MopNumber,Mop2)
    
%matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SM.mat' ];

% load SG struct array
matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SG.mat' ];
if exist(matfile,'file')
load(matfile,'SG');
if size(SG,2) > 0

% divide mop area into 20 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,20,[-100 0]);
    

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



for nsn=1:size(SG,2)
      sn=nsn;
    
      SurveyDatenum=SG(nsn).Datenum;
      %fprintf('Survey %i of %i : %s\n',...
%       nsn,size(SG,2),datestr(SurveyDatenum));
    
    SM(nsn).Mopnum=SG(sn).Mopnum;
    SM(nsn).UTMzone=SG(sn).UTMzone;
    SM(nsn).File=SG(sn).File;
    SM(nsn).FileDatenum=SG(sn).FileDatenum;
    %SM(nsn).Bytes=SG(sn).Bytes;
    SM(nsn).Source=SG(sn).Source;
    SM(nsn).Datenum=SG(sn).Datenum;
    
    SM(nsn).X1D=x1d;
          
% reconstruct x,y grid points within gridded x,y area that 
%  have "no data" NaN's

if ~isempty(SG(sn).X)
xg=SG(sn).X;yg=SG(sn).Y;zg=SG(sn).Z;
cg=SG(sn).Class;
%cg=SG(sn).X*0;
%cmin=accumarray([xidx(:), yidx(:)], cg.',[],@nanmax);
nx=max(xg)-min(xg)+3;ny=max(yg)-min(yg)+3; % temp grid dimensions
tg=nan(ny,nx); % temp grid of Nans
idx=sub2ind([ny nx],yg-min(yg)+2,xg-min(xg)+2); % data point 1d indices
tg(idx)=zg; % add data to temp grid
[X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1);

% now 2d interpolate z values for the Mop transect points
zt=xt*NaN; % initialize transect elevation points 
zt(:) = interp2(X,Y,tg,xt,yt);
  SM(nsn).Z1Dtransect=zt;

% now 2d interpolate z values for all the subtransect points
zst=xst*NaN; % initialize transect elevation points 
zst(:) = interp2(X,Y,tg,xst(:),yst(:));
% get mean, median, std, min-max z(xt) transect values for 51 subtransects.

  SM(nsn).Z1Dmean=nanmean(zst,1);
  SM(nsn).Z1Dmedian=nanmedian(zst,1);
  SM(nsn).Z1Dstd=nanstd(zst,1);
  SM(nsn).Z1Dmin=nanmin(zst,[],1);
  SM(nsn).Z1Dmax=nanmax(zst,[],1);
  
  tg(idx)=cg; % add class data to temp grid
  % now 2d interpolate class values for all the subtransect points
  cst=xst*NaN; % initialize transect class points 
  cst(:) = interp2(X,Y,tg,xst(:),yst(:));
  % find largest class id (hardest substrate)
  SM(nsn).Z1Dclass=nanmax(cst,[],1);

  else
     SM(nsn).Z1Dtransect=[];
     SM(nsn).Z1Dmean=[];
     SM(nsn).Z1Dmedian=[];
     SM(nsn).Z1Dstd=[];
     SM(nsn).Z1Dmin=[];
     SM(nsn).Z1Dmax=[];
     SM(nsn).Z1Dclass=[];
end
  
end

% save survey morpho results in SM file for Mop number

% first sort by date 
if size(SM,2) > 1
T=struct2table(SM); 
sortedT = sortrows(T, 'Datenum');
SM=table2struct(sortedT)';
end

%fprintf('Writing changed SM struct array matfile...\n')
matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SM.mat'];
SMfiles=struct('File',{SM.File},'FileDatenum',num2cell([SM.FileDatenum]));
       eval(['save ' matfile ' SM SMfiles']);
end
end
  else
      fprintf('Skipping creation of an SM file for Mop 1 at the border.\n')
  end
 end

end