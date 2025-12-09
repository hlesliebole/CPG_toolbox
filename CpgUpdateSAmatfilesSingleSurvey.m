function [Survey]=CpgUpdateSAmatfilesSingleSurvey(Survey,SurvIdx)

% Updates all the SA matfiles of mops with data in the 
%  survey file Survey(SurveyIdx).File
%
%  For each mop in Survey(SurveyIdx).NearestMops, it deletes 
%  any existing entry for that survey in the SA 
%  array before adding the new version of the files data. 
%   

CpgDefineMopPath

%  Uses: MopTableUTM.mat, readSurveyFileUTM.m XY2MopNumsV2.m

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');

% initialize individual mop area variables
cmop=[];cx=[];cy=[];cdist=[];

% check to make sure group and drone folders on reefbreak are mounted
%if exist('/volumes/group','dir')  && exist('/volumes/drone','dir')
if exist('/volumes/group','dir')  
    
%----------------------------------------
% read in Survey x,y,z,c (classification) data in UTM coords
%----------------------------------------
 %[xutm,yutm,z,c]=readSurveyFileUTM(Survey(SurvIdx).Source,Survey(SurvIdx).File); 
[xutm,yutm,z,c,utmzone]=readSurveyFileUTM2(Survey(SurvIdx).Source,Survey(SurvIdx).File);
 %fprintf('Min-Max Classification Values: %d %d\n',min(c),max(c))
%----------------------------------------
% reduce to 1m spatial averages  
%----------------------------------------
 
Res=1; % 1m spatial resolution

% round survey x,y to desired resolution
xr=Res*round(xutm/Res); % round to Res meter spatial resolution
yr=Res*round(yutm/Res); % 

% bin and average rounded survey data by placing in unique
%  x,y data array
[ux, ~, xidx] = unique(xr);
[uy, ~, yidx] = unique(yr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray([xidx(:), yidx(:)], 1);  
%array of average of z that fall into each unique x/y combination
zavg = accumarray([xidx(:), yidx(:)], z.')./zcount;
cmode = accumarray([xidx(:), yidx(:)], c.',[], @mode); % most common class 
%tavg = accumarray([xidx(:), yidx(:)], t.')./zcount;
%create a list of the z that fall into each unique x/y combination
%zs = accumarray([xidx(:), yidx(:)], z.', [], @(V) {V}, {});

% reduce arrays to 1d vectors of x,y points with z data 
ii=isnan(zavg(:)) == 0; % 1d indices of valid data
[i,j]=find(isnan(zavg) == 0); % 2d indices of valid data
% final shore box data vectors
xutm=ux(i);yutm=uy(j);
zavg=zavg(ii);
cmode=cmode(ii);
%tavg=tavg(ii);

%fprintf('Min-Max Classification Mode: %d %d\n',min(cmode),max(cmode))

%--------------------------------------------------------
% option to remove spatially averaged data based on less 
%  than N survey points
%
% zcount=zcount(ii); 
% N=3;
%i=find(zcount < N);xutm(i)=[];yutm(i)=[];zavg(i)=[];
%  zcount(i)=[];cmode(i)=[];
%--------------------------------------------------------

%end % end if for npts > 0

fprintf(1,...
    'Reduced to %g , %g x %g meter spatially averaged survey points\n',...
    length(xutm),Res,Res);

%---------------------------------------------------------------- 
% isolate mop index numbers that fall within the survey area bounds
%---------------------------------------------------------------- 

if ~isempty(xutm)
  [nmop]=XY2MopNumsV2(xutm,yutm,Mop);
  umop=unique(sort(nmop));
else
  umop=[];
end

% compare the 1-m avg unique mops to the Survey NearestMops
cidx=find(~ismember([Survey(SurvIdx).NearestMops],umop));
 if numel(cidx) > 0
%      cidx
%      Survey(SurvIdx).NearestMops
%      umop
     fprintf('Mop %i does not have 1m avg data. Updating NearestMops field in Survey index %i.\n');
     Survey(SurvIdx).NearestMops=umop;
     %umop
     %save('SurveyMasterListWithMops.mat','Survey');
 end
 %------------------------------------------------------------------
 % assign survey points to individual mop structural arrays
 %------------------------------------------------------------------
 
 % initialiaze multi SA struct array
 mn=0; % mop counter
 for nn=umop  % sio 2 scarp experiment mops only
     mn=mn+1;
     M(mn).SA=[];
 end 
 
 % add data to multi struct array
 mn=0; % mop counter
 for nn=umop  % sio 2 scarp experiment mops only
     mn=mn+1;
         
% survey counter
 if isempty(M(mn).SA)
    mps=1; % new struct array  
 else   
    mps=size(M(mn).SA,2)+1; 
 end
 

     % survey point indices for this mop number
     midx=find(nmop == nn);
 %    fprintf(' Mop: %i  pts: %i\n',nn,length(midx))
      if ~isempty(midx)

       % pass survey meta info along to mop struct array
         M(mn).SA(mps).Mopnum=nn;
         M(mn).SA(mps).UTMzone='11 S';
         M(mn).SA(mps).File=Survey(SurvIdx).File;
         M(mn).SA(mps).FileDatenum=Survey(SurvIdx).FileDatenum;
         M(mn).SA(mps).Bytes=Survey(SurvIdx).Bytes;
         M(mn).SA(mps).Source=Survey(SurvIdx).Source;
         M(mn).SA(mps).Datenum=Survey(SurvIdx).Datenum;
         M(mn).SA(mps).X=xutm(midx);
         M(mn).SA(mps).Y=yutm(midx);
         M(mn).SA(mps).Z=zavg(midx);
         M(mn).SA(mps).Class=cmode(midx);

     end
 end


fprintf('%s\n','Updating mat files...')
mn=0;
for nn=umop
    mn=mn+1;
    if ~isempty(M(mn).SA)
    matfile=[mpath 'M' num2str(nn,'%5.5i') 'SA.mat'];
      if exist(matfile,'file') % update existing file
        NSA=M(mn).SA;  
        eval(['load ' matfile ' SA']);
        % remove any existing entry for this survey file name
        if size(SA,2) > 0 
         idx=find(strcmp({SA.File},Survey(SurvIdx).File));
        else
            idx=[];
        end
        if numel(idx) > 0
            SA(idx)=[];
        end
        % add new survey entry
        SA=[SA NSA];
        if size(SA,2) > 1
         T=struct2table(SA); % sort by date before saving
         sortedT = sortrows(T, 'Datenum');
         SA=table2struct(sortedT)';
        end
        SAfiles=struct('File',{SA.File},'FileDatenum',...
            num2cell([SA.FileDatenum]));
        eval(['save ' matfile ' SA SAfiles']);
        %fprintf('Mop: %i  NumSurv= %i\n',nn,size(SA,2));
      else % start new file
        SA=M(mn).SA;
%         T=struct2table(SA); % sort by date before saving
%         sortedT = sortrows(T, 'Datenum');
%         SA=table2struct(sortedT)';
        SAfiles=struct('File',{SA.File},'FileDatenum',...
            num2cell([SA.FileDatenum]));
        eval(['save ' matfile ' SA SAfiles']);
        %fprintf('Mop: %i  NumSurv= %i\n',nn,size(SA,2));
      end
    end
end

else
    %fprintf('***** ERROR ******\nEither /volumes/group or /volumes/drone is not mounted.\nStopping.\n')
    fprintf('***** ERROR ******\n/volumes/group is not mounted.\nStopping.\n')
end

end