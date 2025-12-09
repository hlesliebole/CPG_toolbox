% UpdateSAmatfilesMopRange.m
%
%  Build Survey Averaged struct array matfiles for a specified 
%  range of Mops
%
%  Uses: SurveyMasterList.mat, MopTableUTM.mat, readSurveyFileUTM.m 

% set Mop Range to make SA matfiles

%       State Beaches
% Border Field State Park: 1-38
% Silver Strand State Beach: 85-157
% Torrey Pines State Beach: 535-607
% Cardiff State Beach: 664-683
% San Elijo State Beach: 683-706
% Moonlight State Beach: 717-725
% Leucadia State Beach: 740-757
% South Carlsbad State Beach: 762-819
% Carlsbad State Beach: 825-854
% SCARP: 495-636
clear all

Mop1=740;Mop2=854; 

% Load Survey List
load UpdateSurveyList.mat
Survey=Update;

for n=1:size(Update,2)
    Survey(n).MopStart=min([Survey(n).Mops]);
    Survey(n).MopEnd=max([Survey(n).Mops]);
end

% reduce list of survey files to those that have data
%  from Mop1 to Mop2.  Allow for 1 extra mop on the
%  ends of the specified survey mop range in the
%  file naming convention.
%
idx=find([Survey.MopStart]-1 <= Mop2 & [Survey.MopEnd]+1 >= Mop1);
Survey=Survey(idx);

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');

% initialize individual mop area variables
cmop=[];cx=[];cy=[];cdist=[];

% Loop through surveys that include sio scarp mop points

nstep=10; % for updates, just do 1 survey at a time before saving results in the individual
          % mop files.

for loop=1:nstep:size(Survey,2)
  lend=loop+nstep-1;if(lend > size(Survey,2));lend = size(Survey,2);end


 % initial struct array
 mn=0; % mop counter
 for nn=Mop1:Mop2  % sio 2 scarp experiment mops only
     mn=mn+1;
     M(mn).SA=[];
 end 

for SurvNum=loop:lend
 fprintf('Survey number: %i of %i\n',SurvNum,size(Survey,2))
 
 %----------------------------------------
 % read in Survey x,y,z,c (classification) data in UTM coords
 %----------------------------------------
 [xutm,yutm,z,c]=readSurveyFileUTM(Survey(SurvNum).Source,Survey(SurvNum).File); 
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

% %-----------------------------------------------------------------
% 
% %  make interpolated transect line points with approx 1m spatial 
% %  resolution for each Mop within the survey area bounds
% 
% mtx=[];
% mty=[];
% mnum=[];
% extend=100; % extend transect off and onshore, in meters
% 
% 
% if(Survey(SurvNum).MopStart > 1);MopStart=Survey(SurvNum).MopStart-1;else;...
%     MopStart=Survey(SurvNum).MopStart;end
% 
% if(Survey(SurvNum).MopEnd < size(Mop,1));MopEnd=Survey(SurvNum).MopEnd+1;else;...
%     MopEnd=Survey(SurvNum).MopEnd;end
% 
% 
% for MopNum=MopStart:MopEnd
%  % figure out x step size for approx 1m alongshore transect resolution
%  xstep=min([abs(Mop.BackXutm(MopNum)-Mop.OffXutm(MopNum))...
%      abs(Mop.BackYutm(MopNum)-Mop.OffYutm(MopNum))])/...
%         abs(Mop.BackYutm(MopNum)-Mop.OffYutm(MopNum));
% 
%  xi=min([Mop.BackXutm(MopNum) Mop.OffXutm(MopNum)])-extend:xstep:...
%      max([Mop.BackXutm(MopNum) Mop.OffXutm(MopNum)]+extend);
%  yi=interp1([Mop.BackXutm(MopNum) Mop.OffXutm(MopNum)],...
%      [Mop.BackYutm(MopNum) Mop.OffYutm(MopNum)],...
%     xi,'linear','extrap');
% 
%  mtx=[mtx xi]; 
%  mty=[mty yi];
%  mnum=[mnum MopNum*ones(size(xi))];
% end
% 
% %---------------------------------------------------------------------
% 
%  % match survey utm points to the nearest mop transect line
% 
% [dp,lmop]=pdist2([mty',mtx'],[yutm,xutm],'euclidean','smallest',1);

N=XY2MopNums(round(xutm),round(yutm),Mop);
 
 %------------------------------------------------------------------
 % assign survey points to individual mop structural arrays
 %------------------------------------------------------------------
 
 % add to struct array
 mn=0; % mop counter
 for nn=Mop1:Mop2  % defined mop range only
     mn=mn+1;
         
% survey counter
 if isempty(M(mn).SA)
    mps=1; % new struct array  
 else   
    mps=size(M(mn).SA,2)+1; 
 end
 

     % survey point indices for this mop number
     %midx=find(mnum(lmop) == nn);
     midx=find(N == nn);
 %    fprintf(' Mop: %i  pts: %i\n',nn,length(midx))
      if ~isempty(midx)

       % pass survey meta info along to mop struct array
         M(mn).SA(mps).Mopnum=nn;
         M(mn).SA(mps).Datenum=Survey(SurvNum).Datenum;
         M(mn).SA(mps).Source=Survey(SurvNum).Source;
         M(mn).SA(mps).File=Survey(SurvNum).File;
         M(mn).SA(mps).UTMzone='11 S';
         M(mn).SA(mps).X=xutm(midx);
         M(mn).SA(mps).Y=yutm(midx);
         M(mn).SA(mps).Z=zavg(midx);
         M(mn).SA(mps).Class=cmode(midx);

     end
 end
end

fprintf('%s\n','Updating mat files...')
mn=0;
for nn=Mop1:Mop2  % defined mop range only
    mn=mn+1;
    if ~isempty(M(mn).SA)
    matfile=['M' num2str(nn,'%5.5i') 'SA.mat'];
      if exist(matfile,'file') % update existing file
        NSA=M(mn).SA;  
        eval(['load ' matfile ' SA']);
        SA=[SA NSA];
        if(size(SA,2) > 1)
         T=struct2table(SA); % sort by date before saving
         sortedT = sortrows(T, 'Datenum');
         SA=table2struct(sortedT)';
         SA=CombineSurveyDateSourceData(SA);
%             % omit any duplicate survey entries
%             nu=2;
%              while nu <= size(SA,2) 
%                if nu < 3
%                  if strcmp(SA(nu).File,SA(nu-1).File)
%                     SA(nu)=[];
%                  else
%                      nu=nu+1;
%                  end
%         
%                else
%       
%                 if strcmp(SA(nu).File,SA(nu-1).File)...
%                         || strcmp(SA(nu).File,SA(nu-2).File)
%                     SA(nu)=[];
%                 else
%                      nu=nu+1;
%                 end
%                 
%                end
%             end
        end
        eval(['save ' matfile ' SA']);
        fprintf('Mop: %i  NumSurv= %i\n',nn,size(SA,2));
      else % start new file
        SA=M(mn).SA;
        if(size(SA,2) > 1)
         T=struct2table(SA); % sort by date before saving
         sortedT = sortrows(T, 'Datenum');
         SA=table2struct(sortedT)';
         SA=CombineSurveyDateSourceData(SA);
%          % omit any duplicate survey entries
%             nu=2;
%              while nu <= size(SA,2) 
%                if nu < 3
%                  if strcmp(SA(nu).File,SA(nu-1).File)
%                     SA(nu)=[];
%                  else
%                      nu=nu+1;
%                  end
%         
%                else
%       
%                 if strcmp(SA(nu).File,SA(nu-1).File)...
%                         || strcmp(SA(nu).File,SA(nu-2).File)
%                     SA(nu)=[];
%                 else
%                      nu=nu+1;
%                 end
%                 
%                end
%             end
        end
        eval(['save ' matfile ' SA']);
        fprintf('Mop: %i  NumSurv= %i\n',nn,size(SA,2));
      end
    end
end

end