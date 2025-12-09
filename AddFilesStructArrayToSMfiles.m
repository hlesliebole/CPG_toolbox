clearvars

% adds a SAfiles struct array to all the SA matfiles
%mpath='/Volumes/group/MOPS/';
CpgDefineMopPath

load MopTableUTM.mat
load SurveyMasterListWithMops.mat

MopStart=1; % MX border 
MopEnd=11594; % OR border

% Cycle through all existing Mop SA files any struct array
%  entries for file name.
ne=[];
for m=MopStart:MopEnd
    matfile=[mpath 'M' num2str(m,'%5.5i') 'SM.mat'];
    if exist(matfile,'file')
      ll=fprintf('%s\n',matfile);
      load(matfile,'SM');
      if size(SM,2) > 0
       SMfiles=struct('File',{SM.File},'FileDatenum',num2cell([SM.FileDatenum]));
       eval(['save ' matfile ' SM SMfiles']);
      else
       fprintf('** Warning: SM array is empty\n') 
       ne=[ne m];
      end
     end
end
ne