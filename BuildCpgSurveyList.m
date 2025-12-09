%load M00636SA.mat
clear all
CpgSurvey.Mopnum=[];
CpgSurvey.Datenum=[];
CpgSurvey.Source=[];
CpgSurvey.File=[];


addpath ..
fls=dir('../M*SA.mat');

for nn=1:size(fls,1)

%matfile='M00683SA.mat';
matfile=fls(nn).name;
fprintf(' -- %s  --\n',matfile)
%matfile=['M' num2str(nn,'%5.5i') 'SA.mat'];
eval(['load ' matfile ' SA']);

SA=CombineSurveyDateSourceData(SA);

if numel([CpgSurvey.Datenum]) == 0
    
    for n=1:size(SA,2)
          CpgSurvey(n).Mopnum=SA(n).Mopnum;
          CpgSurvey(n).Datenum=SA(n).Datenum;
          CpgSurvey(n).Source=SA(n).Source;
          CpgSurvey(n).File=SA(n).File;
    end
    
else 
    
    n1=size(CpgSurvey,2);

   for n=1:size(SA,2)
%       idx=find([CpgSurvey.Datenum] == SA(n).Datenum & ...
%           strcmpi({CpgSurvey.Source}, SA(n).Source)==1 );
%       if ~isempty(idx)
          CpgSurvey(n1+n).Mopnum=SA(n).Mopnum;
          CpgSurvey(n1+n).Datenum=SA(n).Datenum;
          CpgSurvey(n1+n).Source=SA(n).Source;
          CpgSurvey(n1+n).File=SA(n).File;
 %     end
   end
    
end

end

fprintf('\n Saving CpgSurvey to CpgSurveyMasterList.mat\n')
save CpgSurveyMasterList.mat CpgSurvey
