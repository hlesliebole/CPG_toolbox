clear all
%MakeSurveyMasterList
load SurveyMasterList.mat

%load TestSA.mat
MopNumber=2;

MopPath='/volumes/group/MOPS/';
matfile=[MopPath 'M' num2str(MopNumber,'%5.5i') 'SA.mat'];
load(matfile,'SA');


% first see if there are file name in the SA file that
%  are not found in the Survey struct array (ie. reefbreak1
%  database file names have changed or been deleted for some reason)
idx0=find(~ismember(lower({SA.File}),lower({Survey.File})));

fprintf('-----------------------------------------------------\n')
fprintf('There are %i File names in SA that are not in Survey.\n',numel(idx0))
fprintf('-----------------------------------------------------\n')

% loop through and see if there is another Survey struct file
%  for the same date and source.

for n=idx0
     
    idx2=find( [Survey.Datenum] == SA(n).Datenum & ...
        strcmpi({Survey.Source},SA(n).Source) == 1  &...
         [Survey.Bytes] == SA(n).Bytes );
%     [Survey.MopStart] <= MopNumber &...
%     [Survey.MopEnd] >= MopNumber);
    if ~isempty(idx2)
        if numel(idx2) > 1
            fprintf('**** Multiple Files Same Day: %s \n',datestr(SA(n).Datenum))
            for jj=idx2
                fprintf('%s\n',Survey(jj).File)
            end
        end
        %fprintf('Survey %4i: %s\n',n,Survey(n).File)
        for jj=idx2
            
            fprintf('Old Survey Filename in SA array %4i:\n %s\n',n,SA(n).File)
            fprintf('Date Seems to match Survey      %4i:\n %s\n',jj,Survey(jj).File)
            fprintf('Checking to see if it has the same number of data points...\n')
            [xutm,yutm,z,c]=readSurveyFileUTM(Survey(jj).Source,Survey(jj).File);
            % reduce to 1m spatial averages  
            Res=1; % 1m spatial resolution
            % round survey x,y to desired resolution
            xr=Res*round(xutm/Res);yr=Res*round(yutm/Res); 
            % bin and average rounded survey data by placing in unique x,y data array
            [ux, ~, xidx] = unique(xr);[uy, ~, yidx] = unique(yr);
            %array of counts of the number of points at each unique x/y combination
            zcount = accumarray([xidx(:), yidx(:)], 1);  
            %array of average of z that fall into each unique x/y combination
            zavg = accumarray([xidx(:), yidx(:)], z.')./zcount;
            cmode = accumarray([xidx(:), yidx(:)], c.',[], @mode); % most common class 
            % reduce arrays to 1d vectors of x,y points with z data 
            ii=isnan(zavg(:)) == 0; % 1d indices of valid data
            [i,j]=find(isnan(zavg) == 0); % 2d indices of valid data
            % final shore box data vectors
            Xutm=ux(i);Yutm=uy(j);
            zavg=zavg(ii);
            cmode=cmode(ii);
            % find 1m avg points for this mop number
            Nmop=FindNearestMopTransectsUTM(Xutm,Yutm);
            Xutm=Xutm(Nmop == MopNumber);
            Yutm=Yutm(Nmop == MopNumber);
            zavg=zavg(Nmop == MopNumber);
            cmode=cmode(Nmop == MopNumber);
           
            if numel(SA(n).Z) == numel(zavg)
                % updating file name in SA struct array
                fprintf('--Same number of points. Updating Filename in SA struct array.\n')
                SA(n).File=Survey(jj).File;
            else
                fprintf('**--Different number of points. Adding to end of SA struct array.\n')
                ntot=size(SA,2)+1;
                SA(ntot).Mopnum=MopNumber;
                SA(ntot).Datenum=SA(n).Datenum;
                SA(ntot).File=Survey(jj).File;      
                SA(ntot).X=Xutm;
                SA(ntot).Y=Yutm;
                SA(ntot).Z=zavg;
                SA(ntot).Class=cmode;
            end
            
        end
        
    else
        
        fprintf('****** Survey File Missing ********:\n %s\n',SA(n).File)
        
     end
    
end

% sort by date 
T=struct2table(SA); 
sortedT = sortrows(T, 'Datenum');
SA=table2struct(sortedT)';

%  Now check again if there are file names in the SA file that
%  are not found in the Survey struct array 
idx0=find(~ismember(lower({SA.File}),lower({Survey.File}))); 
% If found delete these entries from the SA struct array
if ~isempty(idx0)
fprintf('-----------------------------------------------------\n')
fprintf('Removing data for %i File names in SA that are not in Survey.\n',numel(idx0))
fprintf('-----------------------------------------------------\n')
SA(idx0)=[]; 
end

fprintf('\n-------------------------------------------------------------\n')
fprintf('-------------------------------------------------------------\n\n')


% find Survey file names not in the SA struct array    
idx=find(~ismember(lower({Survey.File}),lower({SA.File})) &...
    [Survey.MopStart] <= MopNumber &...
    [Survey.MopEnd] >= MopNumber);

fprintf('-----------------------------------------------------\n')
fprintf('There are %i File names in Survey that are not in SA.\n',numel(idx))
fprintf('-----------------------------------------------------\n')

% loop through these surveys and see if the dates and
%   source name actually match a SA entry (eg. file name
%   change or the file was combined with another to form
%   a single survey.

for n=idx
    fprintf('Adding:\n %s\n',Survey(n).File)
            [xutm,yutm,z,c]=readSurveyFileUTM(Survey(n).Source,Survey(n).File);
            % reduce to 1m spatial averages  
            Res=1; % 1m spatial resolution
            % round survey x,y to desired resolution
            xr=Res*round(xutm/Res);yr=Res*round(yutm/Res); 
            % bin and average rounded survey data by placing in unique x,y data array
            [ux, ~, xidx] = unique(xr);[uy, ~, yidx] = unique(yr);
            %array of counts of the number of points at each unique x/y combination
            zcount = accumarray([xidx(:), yidx(:)], 1);  
            %array of average of z that fall into each unique x/y combination
            zavg = accumarray([xidx(:), yidx(:)], z.')./zcount;
            cmode = accumarray([xidx(:), yidx(:)], c.',[], @mode); % most common class 
            % reduce arrays to 1d vectors of x,y points with z data 
            ii=isnan(zavg(:)) == 0; % 1d indices of valid data
            [i,j]=find(isnan(zavg) == 0); % 2d indices of valid data
            % final shore box data vectors
            Xutm=ux(i);Yutm=uy(j);
            zavg=zavg(ii);
            cmode=cmode(ii);
            % find 1m avg points for this mop number
            Nmop=FindNearestMopTransectsUTM(Xutm,Yutm);
            Xutm=Xutm(Nmop == MopNumber);
            Yutm=Yutm(Nmop == MopNumber);
            zavg=zavg(Nmop == MopNumber);
            cmode=cmode(Nmop == MopNumber);
             % add to end of SA array
                ntot=size(SA,2)+1;
                SA(ntot).Mopnum=MopNumber;
                SA(ntot).Datenum=Survey(n).Datenum;
                SA(ntot).File=Survey(n).File;      
                SA(ntot).X=Xutm;
                SA(ntot).Y=Yutm;
                SA(ntot).Z=zavg;
                SA(ntot).Class=cmode;
            
end

% sort by date 
T=struct2table(SA); 
sortedT = sortrows(T, 'Datenum');
SA=table2struct(sortedT)';
