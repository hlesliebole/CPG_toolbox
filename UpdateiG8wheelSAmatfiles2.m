% UpdateSAmatfiles.m
%
%  Updates Survey Averaged struct array matfiles for each Mop area
%
%  Uses: NewSurveySIO2SCARP.mat, MopTableUTM.mat, readSurveyFileUTM.m 

% folder containing CPG Mop mat files

MopPath='/volumes/group/MOPS/';

% Load iG8wheel Survey List

MakeiG8wheelSurveyList

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');

% initialize individual mop area variables
cmop=[];cx=[];cy=[];cdist=[];

% Loop through surveys 

ncs=6; % number of most recent surveys to check if new

for SurvNum=size(Survey,2)-ncs+1:size(Survey,2)
    
 fprintf('Survey number: %g \n',SurvNum)
    Mop1=Survey(SurvNum).MopStart;
    Mop2=Survey(SurvNum).MopEnd;
    
% check if survey has already been added to SA files
% for all the mops in the files mop range
    sdone=1; % set default flag to survey has already been added
    mopnum=Mop1;
    while sdone == 1 && mopnum <= Mop2
        
        matfile=[MopPath 'M' num2str(mopnum,'%5.5i') 'SA.mat'];
        if exist(matfile,'file') % update existing file  
        eval(['load ' matfile ' SA']);
        
        if isfield(SA,'Bytes')
         fprintf(' SA has Bytes field.\n')
        else
         fprintf(' Adding SA Bytes field values...\n') 
          for n=1:size(SA,2)
           sfile=dir(SA(n).File);
           if isempty(sfile)
               SA(n).Bytes=0;
           else
               SA(n).Bytes=sfile.bytes;
           end
          end
          eval(['save ' matfile ' SA']);
        end
      
        
         idx=find(contains({SA.File},Survey(SurvNum).File,'IgnoreCase',true) == 1);
         if isempty(idx) % no survey file match in SA array
             fprintf('Survey Not Found in: %s \n',matfile)
             sdone=0; % set flag to not processed yet
         end
        else
          fprintf(' %s does not exist yet.\n',matfile)
             sdone=0; % set flag to not processed yet  
        end
        
        mopnum=mopnum+1;
    end
    
 if sdone == 1
        fprintf(' Survey has already been included.\n')
 else
    
    % initial struct array
    mn=0; % mop counter
    for nn=Mop1:Mop2  % single file mop range only
       mn=mn+1;
       M(mn).SA=[];
    end 
 
 %----------------------------------------
 % read in Survey x,y,z,c (classification) data in UTM coords
 %----------------------------------------
 [xutm,yutm,z,c]=readSurveyFileUTM(Survey(SurvNum).Source,Survey(SurvNum).File); 
 fprintf('Elevation Range: %6.2f %6.2f\n',min(z),max(z))
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

 %mopidx=Survey(SurvNum).MopStart:Survey(SurvNum).MopEnd;

%-----------------------------------------------------------------

%  make interpolated transect line points with approx 1m spatial 
%  resolution for each Mop within the survey area bounds

mtx=[];
mty=[];
mnum=[];
extend=100; % extend transect off and onshore, in meters


if(Survey(SurvNum).MopStart > 1);MopStart=Survey(SurvNum).MopStart-1;else;...
    MopStart=Survey(SurvNum).MopStart;end

if(Survey(SurvNum).MopEnd < size(Mop,1));MopEnd=Survey(SurvNum).MopEnd+1;else;...
    MopEnd=Survey(SurvNum).MopEnd;end


for MopNum=Mop1:Mop2
 % figure out x step size for approx 1m alongshore transect resolution
 xstep=min([abs(Mop.BackXutm(MopNum)-Mop.OffXutm(MopNum))...
     abs(Mop.BackYutm(MopNum)-Mop.OffYutm(MopNum))])/...
        abs(Mop.BackYutm(MopNum)-Mop.OffYutm(MopNum));

 xi=min([Mop.BackXutm(MopNum) Mop.OffXutm(MopNum)])-extend:xstep:...
     max([Mop.BackXutm(MopNum) Mop.OffXutm(MopNum)]+extend);
 yi=interp1([Mop.BackXutm(MopNum) Mop.OffXutm(MopNum)],...
     [Mop.BackYutm(MopNum) Mop.OffYutm(MopNum)],...
    xi,'linear','extrap');

 mtx=[mtx xi]; 
 mty=[mty yi];
 mnum=[mnum MopNum*ones(size(xi))];
end

%---------------------------------------------------------------------

 % match survey utm points to the nearest mop transect line

[dp,lmop]=pdist2([mty',mtx'],[yutm,xutm],'euclidean','smallest',1);
 
 %------------------------------------------------------------------
 % assign survey points to individual mop structural arrays
 %------------------------------------------------------------------
 
 % add to struct array
 mn=0; % mop counter
 for nn=Mop1:Mop2  % 
     mn=mn+1;
         
% survey counter
 if isempty(M(mn).SA)
    mps=1; % new struct array  
 else   
    mps=size(M(mn).SA,2)+1; 
 end
 

     % survey point indices for this mop number
     midx=find(mnum(lmop) == nn);
 %    fprintf(' Mop: %i  pts: %i\n',nn,length(midx))
      if ~isempty(midx)

       % pass survey meta info along to mop struct array
         M(mn).SA(mps).Mopnum=nn;
         M(mn).SA(mps).Datenum=Survey(SurvNum).Datenum;
         M(mn).SA(mps).Source=Survey(SurvNum).Source;
         M(mn).SA(mps).File=Survey(SurvNum).File;
         M(mn).SA(mps).Bytes=Survey(SurvNum).Bytes;
         M(mn).SA(mps).UTMzone='11 S';
         M(mn).SA(mps).X=xutm(midx);
         M(mn).SA(mps).Y=yutm(midx);
         M(mn).SA(mps).Z=zavg(midx);
         M(mn).SA(mps).Class=cmode(midx);
         
         % correct z datum of 12/6/2021 survey for 2m pole length
         if M(mn).SA(mps).Datenum == datenum(2021,12,6)
             fprintf(1,'****Correcting 12/6/21 survey for pole length\n')
             min(M(mn).SA(mps).Z)
             M(mn).SA(mps).Z=zavg(midx)-4;
             min(M(mn).SA(mps).Z)
         end
      % correct z datum of 1/12/2022 survey for 104.5cm pole length
         if M(mn).SA(mps).Datenum == datenum(2022,1,12)
             fprintf(1,'****Correcting 1/12/22 survey for 104.5cm pole length\n')
%              min(M(mn).SA(mps).Z)
             M(mn).SA(mps).Z=zavg(midx)-0.045;
%              min(M(mn).SA(mps).Z)
         end       

     end
 end


fprintf('%s\n','Updating mat files...')
mn=0;
for nn=Mop1:Mop2 % sio 2 scarp experiment mops only
    mn=mn+1;
    if ~isempty(M(mn).SA)
    matfile=[MopPath 'M' num2str(nn,'%5.5i') 'SA.mat'];
      if exist(matfile,'file') % update existing file
        NSA=M(mn).SA;  
        eval(['load ' matfile ' SA']);
        % skip if survey already exists at this mop
        idx=find(contains({SA.File},Survey(SurvNum).File,'IgnoreCase',true) == 1);
        if isempty(idx) % no survey file match in SA array   
        SA=[SA NSA];
        T=struct2table(SA); % sort by date before saving
        sortedT = sortrows(T, 'Datenum');
        SA=table2struct(sortedT)';
        eval(['save ' matfile ' SA']);
        fprintf('Mop: %i  NumSurv= %i\n',nn,size(SA,2));
        end
      else % start new file
        SA=M(mn).SA;
        T=struct2table(SA); % sort by date before saving
        sortedT = sortrows(T, 'Datenum');
        SA=table2struct(sortedT)';
        eval(['save ' matfile ' SA']);
        fprintf('Mop: %i  NumSurv= %i\n',nn,size(SA,2));
      end
    end
end

end
end
