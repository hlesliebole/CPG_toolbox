clear all
%MakeSurveyMasterList
load SurveyMasterList.mat

load M00582SA.mat
MopNumber=582;

% first see if there are file name in the SA file that
%  are not found in the Survey struct array (ie. reefbreak1
%  database file names have changed or been deleted for some reason)
idx0=find(~ismember(lower({SA.File}),lower({Survey.File})));

% loop through and see if there is another Survey struct file
%  for the same date and source.

for n=idx0
     
    idx2=find( [Survey.Datenum] == SA(n).Datenum & ...
        strcmpi({Survey.Source},SA(n).Source) == 1  &...
    [Survey.MopStart] <= MopNumber &...
    [Survey.MopEnd] >= MopNumber);
    if ~isempty(idx2)
        %fprintf('Survey %4i: %s\n',n,Survey(n).File)
        for jj=idx2
        %fprintf('SA     %4i: %s\n',j,SA(j).File)
        % see if this SA file is also in the Survey struct array
%         idx3=find(ismember(lower({Survey.File}),lower({SA(j).File})));
%         if ~isempty(idx3)
%             %fprintf('Matches Survey index %4i\n',idx3)
%         else
            
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
                fprintf('--Different number of points. Replacing data in SA struct array.\n')
                SA(n).File=Survey(jj).File;      
                SA(n).X=Xutm;
                SA(n).Y=Yutm;
                SA(n).Z=zavg;
                SA(n).Class=cmode;
            end
            
        end
        
    else
        
        fprintf('Survey File Missing?: %s\n',SA(n).File)
        fprintf('Checking to see if it has data for this Mop...\n',Survey(n).File)
    end
    
end


fprintf('-------------------------------------------------------------\n')
fprintf('-------------------------------------------------------------\n')
fprintf('-------------------------------------------------------------\n')


% find Survey file names not in the SA struct array    
    idx=find(~ismember(lower({Survey.File}),lower({SA.File})) &...
    [Survey.MopStart] <= MopNumber &...
    [Survey.MopEnd] >= MopNumber);


% loop through these surveys and see if the dates and
%   source name actually match a SA entry (eg. file name
%   change or the file was combined with another to form
%   a single survey.

for n=idx
    
    idx2=find( [SA.Datenum] == Survey(n).Datenum & ...
        strcmpi({SA.Source},Survey(n).Source) == 1);
    if ~isempty(idx2)
        %fprintf('Survey %4i: %s\n',n,Survey(n).File)
        for j=idx2
        %fprintf('SA     %4i: %s\n',j,SA(j).File)
        % see if this SA file is also in the Survey struct array
        idx3=find(ismember(lower({Survey.File}),lower({SA(j).File})));
        if ~isempty(idx3)
            %fprintf('Matches Survey index %4i\n',idx3)
        else
            fprintf('Survey %4i: %s\n',n,Survey(n).File)
            fprintf('SA     %4i: %s\n',j,SA(j).File)
        end
        end
        
    else
        
        fprintf('New Survey : %s\n',Survey(n).File)
        fprintf('Checking to see if it has data for this Mop...\n',Survey(n).File)
    end
    
end

