% Step2getSurveyFileList.m

% Uses "DataSets" struct array of survey file name conventions definitions 
% to generate the struct array "Survey", a master list containing individual 
% survey file information.

% Survey(n).DataSetNum : Index number for DataSets struct array info
% Survey(n).file : Full filename path
% Survey(n).datenum : Survey date number (convenient for plotting)
% Survey(n).source :  Survey (eg. Gps, UTAir, Truck...)
% Survey(n).type : Parsed piece of filename after the mop range info;
%                  eg. can be used to find types that include "jumbo".
% Survey(n).year  : Survey year (for convenience when writing search code)
% Survey(n).month : Survey month (for convenience when writing search code)
% Survey(n).day : Survey day (for convenience when writing search code)
% Survey(n).MopStart : Starting Mop (for convenience when writing search code)
% Survey(n).MopEnd : Ending Mop (for convenience when writing search code)


% Note:  This script contains some logic to deal with redundancies
% and naming inconsistencies that developed over time with the atv-jumbo
% survey database, eg. past reprocessing of some of these surveys to update 
% the datums etc.  This code may require changes (simplification!if this 
% issue is dealt with by moving/changing filenames in the data base at 
% some point.

% Get survey data set info and filename definitions
Step1defineDataSets

% Settings to limit the Mop range to include in the master list

MopStart=495; % start south end LJ shores 
MopEnd=636; % end Del Mar inlet
ns=0; % survey counter

nns=size(DataSets,2); % number of file data sets

clear Survey

for nn=[1 2 4] % limit to SIO Gps and Truck Lidar (data sets 1-3)
%for nn=1:nns % loop through data sets
    
fprintf(1,'Getting %s surveys with wildcard path:\n',DataSets(nn).name)
fprintf(1,'%s\n',DataSets(nn).filepath)
filesd=dir(DataSets(nn).filepath);
ndirs=size(filesd,1);
fprintf('%g survey directories found. Checking for Shore Box overlap...\n',ndirs) 

% loop through directory names, parse out survey date and mop range, and
%  test for overlap with the mop range settings

for n=1:ndirs
    
    %info=regexprep([filesd(n).folder '/' filesd(n).name],DataSets(nn).filehead, '');
    info=regexprep(fullfile(filesd(n).folder,filesd(n).name),...
        regexprep(DataSets(nn).filehead,'\\','\\\'), '');
    fprintf('%s\n',info)
    info=regexprep(info,'_D_', ' ');
    info=regexprep(info,'_', ' ');
    info=regexp(info,'\S*','match');% turns info into cell array of strings
    
    if(size(info,2) > 3)  % check to make sure directory name is legit
        
    SurveyDatenum=datenum(info{1},'yyyymmdd'); % survey date as datenum
    SurveyStart=str2num(info{2}); % survey mop start as number
    SurveyEnd=str2num(info{3}); % survey mop end as number
    SurveyType=strcat(info{4:end});
    %fprintf('%g %g %g %g %g\n',SurveyDatenum,SurveyStart,SurveyEnd,...
    %    MopStart,MopEnd)
    
    % check if survey and mop range overlap
    
    if( (SurveyStart < MopEnd) && SurveyEnd > MopStart) 
        
      %fprintf('Including %s\n',filesd(n).name)
      % if there is overlap, save survey data file name and info
      sfile=dir([filesd(n).folder '/' filesd(n).name]);
      
      fprintf('** Matched ** %g %g %s\n',size(sfile),...
          fullfile(sfile(1).folder,sfile(1).name))
      %[sfile(1).folder '/' sfile(1).name])
      
      if(size(sfile,1) == 1) % check if found a file
          
          % If this is Gps data, then check if the Gps survey date had
          % been previously found (as a .adj2011 file) if so, overwrite
          % that entry with the new (.navd88 file) info.  Otherwise
          % advance the counter and add the new survey to the list
       if(strcmpi(DataSets(nn).name,'gps') && ns > 0)
            igps=find(strcmp({Survey.source}, 'Gps') == 1 );
            i=find(vertcat(Survey(igps).datenum) == SurveyDatenum);
            if(~isempty(i))
            Survey(igps(i)).file=fullfile(sfile(1).folder,sfile(1).name);
            fprintf('** Replaced Gps filename ** %g %s\n',i,Survey(i).file)
            else
            ns=ns+1;     
        % if multiple datafiles in same directory, use first one
        %fprintf('   survey file: %s\n',sfile(1).name)
            Survey(ns).DataSetNum=nn;
            Survey(ns).file=fullfile(sfile(1).folder,sfile(1).name);
            Survey(ns).datenum=SurveyDatenum;
            Survey(ns).source=DataSets(nn).name;
            Survey(ns).type=SurveyType;
            Survey(ns).year=str2num(datestr(SurveyDatenum,'yyyy'));
            Survey(ns).month=str2num(datestr(SurveyDatenum,'mm'));
            Survey(ns).day=str2num(datestr(SurveyDatenum,'dd'));
            Survey(ns).MopStart=SurveyStart;
            Survey(ns).MopEnd=SurveyEnd;
         
          
        fprintf('** Matched ** %g %s\n',ns,Survey(ns).file)
            end
       else
              ns=ns+1;     
        % if multiple datafiles in same directory, use first one
        %fprintf('   survey file: %s\n',sfile(1).name)
            Survey(ns).DataSetNum=nn;
            Survey(ns).file=fullfile(sfile(1).folder,sfile(1).name);
            Survey(ns).datenum=SurveyDatenum;
            Survey(ns).source=DataSets(nn).name;
            Survey(ns).type=SurveyType;
            Survey(ns).year=str2num(datestr(SurveyDatenum,'yyyy'));
            Survey(ns).month=str2num(datestr(SurveyDatenum,'mm'));
            Survey(ns).day=str2num(datestr(SurveyDatenum,'dd'));
            Survey(ns).MopStart=SurveyStart;
            Survey(ns).MopEnd=SurveyEnd;

          
        fprintf('** Matched ** %g %s\n',ns,Survey(ns).file)  
       end
      end
          
    end
    
    else % end bad directory name check
       %fprintf('   *** Bad Survey Directory Name: %s\n',filesd(n).name) 
    end  
end

fprintf('%g surveys found between Mops %g - %g\n',ns,MopStart,MopEnd)
end

Survey2=Survey; % save new survey list as Survey2

% load original survey file Survey struct array
load SurveySIO2SCARP.mat

% reduce Survey2 struct array to new surveys only
for n=size(Survey2,2):-1:1
    idx=find( [Survey.datenum] == Survey2(n).datenum & ...
        strcmpi({Survey.source}, Survey2(n).source )==1);
    if isempty(idx)
        fprintf('%s %s\n',datestr(Survey2(n).datenum),Survey2(n).source)
    else
        Survey2(n)=[];
    end
end

Survey=Survey2; % return to original Survey struct array name

fprintf('Saving new Survey struct array to NewSurveySIO2SCARP.mat\n')
save NewSurveySIO2SCARP.mat Survey