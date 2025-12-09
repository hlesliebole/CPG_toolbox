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
ne=[];
for m=MopStart:MopEnd
    matfile=[mpath 'M' num2str(m,'%5.5i') 'SG.mat'];
    if exist(matfile,'file')
      ll=fprintf('%s\n',matfile);
      load(matfile,'SG');
      if size(SG,2) > 0
       SGfiles=struct('File',{SG.File},'FileDatenum',num2cell([SG.FileDatenum]));
       eval(['save ' matfile ' SG SGfiles']);
      else
       fprintf('** Warning: SG array is empty\n') 
       ne=[ne m];
      end
     end
end
ne