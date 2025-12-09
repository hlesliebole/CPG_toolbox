clear all
%  ReProcessSingleSurvey.m
%mpath='/Volumes/group/MOPS/';
DefineMopPath

load MopTableUTM.mat
load SurveyMasterListWithMops.mat

MopStart=659;%1; % MX border 
MopEnd=11594; % OR border

% Cycle through all existing Mop SA files any struct array
%  entries for file name.

for m=MopStart:MopEnd
    matfile=[mpath 'M' num2str(m,'%5.5i') 'SG.mat'];
    if exist(matfile,'file')
      ll=fprintf('%s\n',matfile);
      load(matfile,'SG');
      if size(SG,2) > 0
%       if isfield(SG,'Bytes')
%        SG=rmfield(SG,'Bytes');% remove the bytes field
%       end
      
      if ~isfield(SG,'FileDatenum')
       SG(1).FileDatenum=[];% initialize the FileDatenum field
      end
      if ~isfield(SG,'Bytes')
       SG(1).Bytes=[];% initialize the FileDatenum field
      end
      
      % add file modified datenums 
      for n=1:size(Survey,2)
       idx=find(strcmpi({SG.File},Survey(n).File) == 1);
       if ~isempty(idx)
          SG(idx).FileDatenum=Survey(n).FileDatenum;
          SG(idx).Bytes=Survey(n).Bytes;
       end
      end
      
      % reorder fields to match the Survey2SA.m database rebuild order
      forder={'Mopnum','UTMzone','File','FileDatenum','Bytes',...
      'Source','Datenum','X','Y','Z','Class'};
  
      SG=orderfields(SG,forder);
          
      save(matfile,'SG');
       
      end
    end
end
