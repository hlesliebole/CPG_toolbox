
% check to make sure group and drone folders on reefbreak are mounted
%if exist('/volumes/group','dir')  && exist('/volumes/drone','dir')
if exist('/volumes/group','dir')  

%--------------
%  Step 1
%---------------------

% Makes "DataSets" struct array with key data set information needed
%  to fine the data files.

DefineCliffDataSets;

%-------------------------------------------------
% Step2
%-------------------------------------------------

% Uses "DataSets" struct array of survey file name conventionsand  definitions 
% to generate the struct array "Survey", a master list containing all
% individual survey file information.

% CliffSurvey(n).DataSetNum : Index number for DataSets struct array info
% CliffSurvey(n).file : Full filename path
% CliffSurvey(n).datenum : Survey date number (convenient for plotting)
% CliffSurvey(n).source :  Survey (eg. Gps, UTAir, Truck...)
% CliffSurvey(n).type : Parsed piece of filename after the mop range info;
%                  eg. can be used to find types that include "jumbo".
% CliffSurvey(n).year  : Survey year (for convenience when writing search code)
% CliffSurvey(n).month : Survey month (for convenience when writing search code)
% CliffSurvey(n).day : Survey day (for convenience when writing search code)
% CliffSurvey(n).MopStart : Starting Mop (for convenience when writing search code)
% CliffSurvey(n).MopEnd : Ending Mop (for convenience when writing search code)


% Note:  This script contains some logic to deal with redundancies
% and naming inconsistencies that developed over time with the atv-jumbo
% survey database, eg. past reprocessing of some of these surveys to update 
% the datums etc.  This code may require changes (simplification! if this 
% issue is dealt with by moving/changing filenames in the data base at 
% some point.

% Settings to limit the Mop range to include in the master list

% MopStart=495; % start south end LJ shores 
% MopEnd=636; % end Del Mar inlet
MopStart=1; % MX border 
MopEnd=11594; % OR border
MopStart=654; % 
MopEnd=654; % 
ns=0; % survey counter

nns=size(DataSets,2); % number of file data sets

if exist('CliffSurveyMasterList.mat','file')
    load CliffSurveyMasterList.mat
    OldCliffSurvey=CliffSurvey;
else
    OldCliffSurvey=[];
end

clear CliffSurvey

for nn=1:nns % loop through data sets
%    for nn=1:2 % loop through data sets
    
fprintf(1,'Getting %s surveys with wildcard path: %s\n',DataSets(nn).name,DataSets(nn).filepath)
%fprintf(1,'%s\n',DataSets(nn).filepath)
filesd=dir(DataSets(nn).filepath);
ndirs=size(filesd,1);
%fprintf('%g survey directories found. Checking for Filename naming convention or date errors...\n\n',ndirs) 

% loop through directory names, parse out survey date and mop range, and
%  test for overlap with the mop range settings

for n=1:ndirs
    
    info=regexprep(fullfile(filesd(n).folder,filesd(n).name),...
        regexprep(DataSets(nn).filehead,'\\','\\\'), '');
 %%%%%   %fprintf('%s\n',info)
    info=regexprep(info,'_D_', ' ');
    info=regexprep(info,'_', ' ');
    info=regexp(info,'\S*','match');% turns info into cell array of strings
    
    % parse second date from file name
    info2=filesd(n).name;
    info2=regexprep(info2,'filtered_clean','');
    info2=regexprep(info2,'\.',' ');
    info2=regexp(info2,'\S*','match');% turns info into cell array of strings
    FileDatenum=datenum(info2{1},'yyyymmdd'); % survey date as datenum
    
    if(size(info,2) > 3)  % check to make sure directory name is legit
    
    SurveyDatenum=datenum(info{1},'yyyymmdd'); % survey date as datenum
    SurveyStart=str2num(info{2}); % survey mop start as number
    SurveyEnd=str2num(info{3}); % survey mop end as number
    SurveyType=strcat(info{4:end});
    
%     if(SurveyDatenum ~= FileDatenum)
%         fprintf('---------\n')
%         fprintf('WARNING: Survey Filename convention error or Date Mismatch with Folder name\n')
%         fprintf('%s\n',filesd(n).folder)
%         fprintf('%s\n',filesd(n).name)
%         fprintf('---------\n')
%     end
        
    %fprintf('%g %g %g %g %g\n',SurveyDatenum,SurveyStart,SurveyEnd,...
    %    MopStart,MopEnd)
    
    % check if survey and mop range overlap
    
    if( ~isempty(SurveyStart) && SurveyStart < MopEnd && SurveyEnd > MopStart) 
        
      %fprintf('Including %s\n',filesd(n).name)
      % if there is overlap, save survey data file name and info
      sfile=dir([filesd(n).folder '/' filesd(n).name]);
      
 %%%%%%%     fprintf('** Matched ** %g %g %s\n',size(sfile),...
 %%%%%%%         fullfile(sfile(1).folder,sfile(1).name))
      %[sfile(1).folder '/' sfile(1).name])
      
      if(size(sfile,1) == 1) % check if found a file
          
          % If this is Gps data, then check if the Gps survey date had
          % been previously found (as a .adj2011 file) if so, overwrite
          % that entry with the new (.navd88 file) info.  Otherwise
          % advance the counter and add the new survey to the list
       if(strcmpi(DataSets(nn).name,'gps') && ns > 0)
            igps=find(strcmp({CliffSurvey.Source}, 'Gps') == 1 );
            i=find(vertcat(CliffSurvey(igps).Datenum) == SurveyDatenum & ...
                vertcat(CliffSurvey(igps).MopStart) == SurveyStart & ...
                vertcat(CliffSurvey(igps).MopEnd) == SurveyEnd);
            if(~isempty(i))
            CliffSurvey(igps(i)).File=fullfile(sfile(1).folder,sfile(1).name);
  %%%%%%          fprintf('** Replaced Gps filename ** %g %s\n',i,CliffSurvey(i).File)
            else
            ns=ns+1;     
        % if multiple datafiles in same directory, use first one
        %fprintf('   survey file: %s\n',sfile(1).name)
            CliffSurvey(ns).DataSetNum=nn;
            CliffSurvey(ns).File=fullfile(sfile(1).folder,sfile(1).name);
            CliffSurvey(ns).FileDatenum=datenum(sfile(1).date);
            CliffSurvey(ns).Bytes=sfile(1).bytes;
            CliffSurvey(ns).Datenum=SurveyDatenum;
            CliffSurvey(ns).Source=DataSets(nn).name;
            CliffSurvey(ns).Type=SurveyType;
            CliffSurvey(ns).Year=str2num(datestr(SurveyDatenum,'yyyy'));
            CliffSurvey(ns).Month=str2num(datestr(SurveyDatenum,'mm'));
            CliffSurvey(ns).Day=str2num(datestr(SurveyDatenum,'dd'));
            CliffSurvey(ns).MopStart=SurveyStart;
            CliffSurvey(ns).MopEnd=SurveyEnd;
         
          
   %%%%%%%%%  fprintf('** Matched ** %g %s\n',ns,CliffSurvey(ns).File)
            end
       else
              ns=ns+1;     
        % if multiple datafiles in same directory, use first one
        %fprintf('   survey file: %s\n',sfile(1).name)
            CliffSurvey(ns).DataSetNum=nn;
            CliffSurvey(ns).File=fullfile(sfile(1).folder,sfile(1).name);
            CliffSurvey(ns).FileDatenum=datenum(sfile(1).date);
            CliffSurvey(ns).Bytes=sfile(1).bytes;
            CliffSurvey(ns).Datenum=SurveyDatenum;
            CliffSurvey(ns).Source=DataSets(nn).name;
            CliffSurvey(ns).Type=SurveyType;
            CliffSurvey(ns).Year=str2num(datestr(SurveyDatenum,'yyyy'));
            CliffSurvey(ns).Month=str2num(datestr(SurveyDatenum,'mm'));
            CliffSurvey(ns).Day=str2num(datestr(SurveyDatenum,'dd'));
            CliffSurvey(ns).MopStart=SurveyStart;
            CliffSurvey(ns).MopEnd=SurveyEnd;

          
 %%%%%%       fprintf('** Matched ** %g %s\n',ns,CliffSurvey(ns).File)  
       end
      end
          
    end
    
    else % end bad directory name check
       %fprintf('   *** Bad Survey Directory Name: %s\n',filesd(n).name) 
    end  
end

%%%%%fprintf('%g surveys found between Mops %g - %g\n',ns,MopStart,MopEnd)
end

fprintf('-----------------------------------------------------\n')
fprintf('%i survey files found.\n',size(Survey,2))


% sort surveys by date
T=struct2table(CliffSurvey);
sortedT = sortrows(T, 'Datenum');
CliffSurvey=table2struct(sortedT)';
clear T sortedT;

%------------------------------------
% old/new survey file list comparisons
%------------------------------------

if ~isempty(OldCliffSurvey)
    
% find any survey files with a changed modification date 
nmod=[];
for n=1:size(CliffSurvey,2)
    idx=find(strcmpi({OldCliffSurvey.File},CliffSurvey(n).File));
    if ~isempty(idx)
        if OldCliffSurvey(idx).FileDatenum ~= CliffSurvey(n).FileDatenum
            nmod=[nmod n];
        end
    end
end
fprintf('%i files have been modified since the previous master list was made.\n',numel(nmod))
for n=nmod
    fprintf('  %s\n',CliffSurvey(n).File)
end

% find filenames in the OldSurvey struct array that
%  are not found in the new Survey struct array 
idx1=find(~ismember(lower({OldCliffSurvey.File}),lower({CliffSurvey.File})));


fprintf('There are %i filenames in the previous master list that no longer exist.\n',numel(idx1))
%fprintf('-----------------------------------------------------\n')

if numel(idx1) > 0
    for n=idx1
     fprintf('%s\n',OldCliffSurvey(n).File)
    end
end

% find filenames in the OldSurvey struct array that
%  are not found in the new Survey struct array 
idx2=find(~ismember(lower({CliffSurvey.File}),lower({OldCliffSurvey.File})));

%fprintf('-----------------------------------------------------\n')
fprintf('There are %i new filenames that were not in the previous master list.\n',numel(idx2))
fprintf('-----------------------------------------------------\n')

if numel(idx2) > 0
    for n=idx2
     fprintf('%s\n',CliffSurvey(n).File)
    end
end


end
%%%%%%fprintf('Saving Survey struct array to SurveyMasterList.mat\n')
save CliffSurveyMasterList.mat CliffSurvey

else
    %fprintf('***** ERROR ******\nEither /volumes/group or /volumes/drone is not mounted.\nStopping.\n')
    fprintf('***** ERROR ******\nEither /volumes/group is not mounted.\nStopping.\n')
end