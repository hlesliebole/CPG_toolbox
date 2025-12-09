clear all
%  ReProcessSingleSurvey.m
%mpath='/Volumes/group/MOPS/';
DefineMopPath

load MopTableUTM.mat
load SurveyMasterListWithMops.mat

MopStart=1; % MX border 
MopEnd=11594; % OR border

% Cycle through all existing Mop SA files any struct array
%  entries for file name.

for m=MopStart:MopEnd
    matfile=[mpath 'M' num2str(m,'%5.5i') 'SA.mat'];
    if exist(matfile,'file')
      ll=fprintf('%s\n',matfile);
      load(matfile,'SA');
     
%       if isfield(SA,'Bytes')
%        SA=rmfield(SA,'Bytes');% remove the bytes field
%       end
      
      if ~isfield(SA,'FileDatenum')
       SA(1).FileDatenum=[];% initialize the FileDatenum field
      end
      
      if ~isfield(SA,'Bytes')
       SA(1).Bytes=[];% initialize the FileDatenum field
      end
      
      
      % add file modified datenums and bytes fields
      for n=1:size(Survey,2)
       idx=find(strcmpi({SA.File},Survey(n).File) == 1);% & ismember(m,[Survey(n).NearestMops]));
       if ~isempty(idx)
          SA(idx).FileDatenum=Survey(n).FileDatenum;
          SA(idx).Bytes=Survey(n).Bytes;
%        else
%            fprintf('Could not find a Survey name match for:\n%s\n',Survey(n).File)
       end
       
      end
      
      % reorder fields to match the Survey2SA.m database rebuild order
      forder={'Mopnum','UTMzone','File','FileDatenum',...
      'Bytes','Source','Datenum','X','Y','Z','Class'};
  
      SA=orderfields(SA,forder);
          
      save(matfile,'SA');
       
    end
end
