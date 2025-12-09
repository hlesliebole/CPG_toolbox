clearvars

% adds a SAfiles struct array to all the SA matfiles
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
      if size(SA,2) > 0
       SAfiles=struct('File',{SA.File},'FileDatenum',num2cell([SA.FileDatenum]));
       eval(['save ' matfile ' SA SAfiles']);
      else
       fprintf('** Warning: SA array is empty\n') 
       SA=[];
       SAfiles=[];
       eval(['save ' matfile ' SA SAfiles']);
      end
     end
end
