clear all
DefineMopPath

%Mop1=550;Mop2=555;% NB
Mop1=569;Mop2=596;% ruby2d

%       State Beaches
%Mop1=26;Mop2=38;% Border Field State Park: 1-38
% Silver Strand State Beach: 85-157
% Mop1=535;Mop2=607;% Torrey Pines State Beach: 535-607
% Cardiff State Beach: 664-683
% San Elijo State Beach: 683-706
% Moonlight State Beach: 717-725
% Leucadia State Beach: 740-757
% South Carlsbad State Beach: 762-819
% Carlsbad State Beach: 825-854

% load struct array list of moved survey data files if it exists
if exist('MovedSurveyFiles.mat','file')
    load MovedSurveyFiles.mat
end

%MakeSurveyMasterList
load SurveyMasterList.mat

%load TestSA.mat
for MopNumber=Mop1:Mop2
    
%fprintf('-----------------------------------------------------\n')
fprintf('-----------------------------------------------------\n')
fprintf('Updating Mop %i SA matfile.\n',MopNumber)
fprintf('-----------------------------------------------------\n')
%fprintf('-----------------------------------------------------\n')


MopPath='/volumes/group/MOPS/';
matfile=[MopPath 'M' num2str(MopNumber,'%5.5i') 'SA.mat'];
load(matfile,'SA');

cflag=0; % struct array change flag set to no change

% remove any duplicate surveys based on the filename
 [~, idx] = unique({SA.File}, 'stable');
 if numel(idx) < size(SA,2)
     fprintf(' Removing %i Duplicate Survey File Entries from SA...\n',...
         size(SA,2)-numel(idx)) 
      ddx=find(~ismember(1:size(SA,2),idx));
      for n=ddx
          fprintf('%s\n',SA(n).File)
      end
      SA=SA(idx);
      cflag=1; % struct array change flag
 else
     fprintf('No duplicate points.\n')
 end

% check if SA file has Bytes field, if not loop through files
%  and add Bytes field
if isfield(SA,'Bytes')
    fprintf(' SA has Bytes field.\n')
else
   fprintf(' Adding SA Bytes field values...\n') 
   cflag=1; % struct array change flag
   for n=1:size(SA,2)
           sfile=dir(SA(n).File);
           if isempty(sfile)
               SA(n).Bytes=0;
           else
               SA(n).Bytes=sfile.bytes;
           end
   end
end

fprintf('Checking for file name changes...\n')
% first loop through SA file names and update any previously
%  found file name or file location changes
if exist('Moved','var')
    %fprintf('Comparing SA to moved file list...\n')
    nmoved=1;
    mpass=0;
   while nmoved > 0 
    nmoved=0;
    mpass=mpass+1;
    for n=1:size(SA,2)
      idxm=find(strcmpi({Moved.OldName},SA(n).File) == 1);
      if ~isempty(idxm)
       nmoved=nmoved+1;
       fprintf('Renaming %s\n',SA(n).File)
       fprintf('      to %s\n',Moved(idxm(end)).NewName)
       SA(n).File=Moved(idxm(end)).NewName;
       SA(n).Bytes=Moved(idxm(end)).Bytes;
       cflag=1; % struct array change flag
      end
    end
    fprintf('Moved %i filenames on pass %i.\n',nmoved,mpass)
   end % end while
end

% second, see if there are file name in the SA file that
%  are not found in the Survey struct array (ie. reefbreak1
%  database file names have changed or been deleted for some reason)
idx0=find(~ismember(lower({SA.File}),lower({Survey.File})));

%fprintf('-----------------------------------------------------\n')
fprintf('There are %i File names in SA that are not in Survey.\n',numel(idx0))
%fprintf('-----------------------------------------------------\n')

% loop through and see if there is another Survey struct file
%  for the same date and source.

for n=idx0
    
    idx2=find( [Survey.Datenum] == SA(n).Datenum & ...
        strcmpi({Survey.Source},SA(n).Source) == 1  &...
        contains({Survey.File},'cart') == 0  &...
    [Survey.MopStart] <= MopNumber &...
    [Survey.MopEnd] >= MopNumber);
    if ~isempty(idx2)
        if numel(idx2) > 1
            fprintf('**** Multiple Files Same Day & Source: %s \n',datestr(SA(n).Datenum))
            for jj=idx2
                fprintf('%s\n',Survey(jj).File)
            end
        end
        %fprintf('Survey %4i: %s\n',n,Survey(n).File)
        for jj=idx2
            
            fprintf('Old Survey Filename in SA array %4i:\n %s\n',n,SA(n).File)
            fprintf('Date Seems to match Survey      %4i:\n %s\n',jj,Survey(jj).File)         
            fprintf('Checking to see if it has the same number of data points...\n')
            [xutm,yutm,z,c,utmzone]=readSurveyFileUTM2(Survey(jj).Source,Survey(jj).File);
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
                fprintf('--Same number of points. Updating Filename in SA struct array,\n')
                fprintf('  and adding to moved file struct array list.\n');
                if exist('Moved','var');nm=size(Moved,2);else;nm=0;end
                Moved(nm+1).OldName=SA(n).File;
                Moved(nm+1).NewName=Survey(jj).File;
                Moved(nm+1).Bytes=Survey(jj).Bytes;
                SA(n).File=Survey(jj).File;
                SA(n).FileDatenum=Survey(jj).FileDatenum;
                SA(n).Bytes=Survey(jj).Bytes;
                cflag=1; % struct array change flag
            else
                fprintf('**--Different number of points. Adding to end of SA struct array.\n')
                ntot=size(SA,2)+1;
                SA(ntot).Mopnum=MopNumber;
                SA(ntot).UTMzone='11 S';
                SA(ntot).File=Survey(jj).File;
                SA(n).FileDatenum=Survey(jj).FileDatenum;
                SA(ntot).Source=Survey(jj).Source;
                SA(ntot).Datenum=SA(n).Datenum;
                SA(ntot).Bytes=Survey(jj).Bytes;    
                SA(ntot).X=Xutm;
                SA(ntot).Y=Yutm;
                SA(ntot).Z=zavg;
                SA(ntot).Class=cmode;
                cflag=1; % struct array change flag
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
%fprintf('-----------------------------------------------------\n')
fprintf('Removing data for %i File names in SA that are not in Survey.\n',numel(idx0))
%fprintf('-----------------------------------------------------\n')
SA(idx0)=[]; 
cflag=1; % struct array change flag
end

%fprintf('\n-------------------------------------------------------------\n')
%fprintf('-------------------------------------------------------------\n\n')


% find Survey file names not in the SA struct array    
idx=find(~ismember(lower({Survey.File}),lower({SA.File})) &...
    [Survey.MopStart] <= MopNumber &...
    [Survey.MopEnd] >= MopNumber);

%fprintf('-----------------------------------------------------\n')
fprintf('There are %i File names in Survey that are not in SA.\n',numel(idx))
%fprintf('-----------------------------------------------------\n')

% loop through these surveys and see if the dates and
%   source name actually match a SA entry (eg. file name
%   change or the file was combined with another to form
%   a single survey.
nadd=0;
for n=idx
    nadd=nadd+1;
    fprintf('Adding:\n %i of %i %s\n',nadd,numel(idx),Survey(n).File)
            [xutm,yutm,z,c,utmzone]=readSurveyFileUTM2(Survey(n).Source,Survey(n).File);
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
                SA(ntot).UTMzone='11 S'; 
                SA(ntot).File=Survey(n).File; 
                SA(ntot).FileDatenum=Survey(n).FileDatenum;
                SA(ntot).Source=Survey(n).Source;
                SA(ntot).Datenum=Survey(n).Datenum;
                SA(ntot).Bytes=Survey(n).Bytes; 
                SA(ntot).X=Xutm;
                SA(ntot).Y=Yutm;
                SA(ntot).Z=zavg;
                SA(ntot).Class=cmode;
                cflag=1; % struct array change flag
            
end

% output changed struct array if changes were made

if cflag == 1
 
    % remove any duplicate surveys based on the filename
 [~, idx] = unique({SA.File}, 'stable');
 if numel(idx) < size(SA,2)
     fprintf(' Removing %i Duplicate Survey File Entries from SA...\n',...
         size(SA,2)-numel(idx)) 
      ddx=find(~ismember(1:size(SA,2),idx));
      for n=ddx
          fprintf('%s\n',SA(n).File)
      end
      SA=SA(idx);
      cflag=1; % struct array change flag
 else
     fprintf('No duplicate points.\n')
 end
 
 fprintf('Saving updated SA struct array in matfile...\n')   
 % sort by date 
 T=struct2table(SA); 
 sortedT = sortrows(T, 'Datenum');
 SA=table2struct(sortedT)';

 % flag known bad points with SA.Class=-999
    [SA]=SurveySpecialQC(SA); 
 save(matfile,'SA');

 if exist('Moved','var')
     save MovedSurveyFiles.mat Moved
 end
else
    fprintf('No changes to SA struct array.\n')   
end

end


% % find dates with more than 1 survey file
% [~, idx] = unique([SA.Datenum], 'stable');
% DupDates=find(~ismember([1:size(SA,2)],idx));
% [~, idx2] = unique(DupDates, 'stable');
% 
% for n=unique(DupDates)
%     fprintf('%s \n',datestr(SA(n).Datenum)) 
%     idp=find([SA.Datenum] == SA(n).Datenum);
%     for i=idp
%       fprintf('%i  : %s\n',numel(SA(i).X),SA(i).File)  
%     end
% end


% %fix 594 missing struct array info
% load /volumes/group/mops/M00594SA.mat
% for n=1:size(SA,2)
%     idx=find(strcmpi({Survey.File},SA(n).File) == 1);
%     SA(n).Source=Survey(idx).Source;
%     SA(n).Bytes=Survey(idx).Bytes;
%     SA(n).UTMzone='11 S';  
% end
% save /volumes/group/mops/M00594SA.mat SA