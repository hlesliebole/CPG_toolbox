clear all
load SurveyMasterList.mat
load CpgSurveyMasterList.mat
nu=0;

for n=1:size(Survey,2)
    fprintf('%s\n',Survey(n).File)
    mlist=[];
    for j=Survey(n).MopStart:Survey(n).MopEnd
       idx=find([CpgSurvey.Mopnum] == j & ...
           [CpgSurvey.Datenum] == Survey(n).Datenum & ...
           strcmpi({CpgSurvey.Source},Survey(n).Source) == 1);
       if isempty(idx)
           mlist=[mlist j];
       end
    end
    if numel(mlist) > 0
     fprintf('%s\n',num2str(mlist))
     nu=nu+1;
     Update(nu).File=Survey(n).File;
     Update(nu).Source=Survey(n).Source;
     Update(nu).Datenum=Survey(n).Datenum;
     Update(nu).Mops=mlist;
    end
    
end

save UpdateSurveyList.mat Update



