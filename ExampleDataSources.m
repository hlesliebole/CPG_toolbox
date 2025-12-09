
MopFile=['M00618SM.mat'];
load(MopFile)
fprintf('\n %s contains %s surveys.\n\n',MopFile,num2str(size(SM,2)))

SurveySources=unique({SM.Source},'stable'); % unique survey sources
fprintf('Survey Sources: %s \n\n',strjoin(SurveySources,', '))

% Example to get indexes of all Gps surveys
GpsIndexes=find(strcmpi({SM.Source}, 'gps')==1);
% Example to get indexes of Gps Jumbo surveys
JumboIndexes=find(contains({SM.File}, 'jumbo','IgnoreCase',true)==1); 

% Make a table of Survey Source Info
fprintf('%10s %15s %15s \n','Source','# of Surveys','Last Survey')
fprintf('%10s %15s %15s \n','------','------------','-----------')
for N=1:length(SurveySources)
    LastSurvey=find(strcmpi({SM.Source}, SurveySources{N})==1,1,'last');
    if strcmpi(SurveySources{N},'gps')
        LastJumbo=find( contains({SM.File}, 'jumbo','IgnoreCase',true) == 1,1,'last');
        fprintf('%10s %10i %20s  Jumbos: %10i %20s\n',SurveySources{N},...
        length(find(strcmpi({SM.Source}, SurveySources{N})==1)),...
        datestr(SM(LastSurvey).Datenum),...
        length(find(contains({SM.File}, 'jumbo',...
          'IgnoreCase',true)==1)),...
          datestr(SM(LastJumbo).Datenum))
    else
        fprintf('%10s %10i %20s\n',SurveySources{N},...
        length(find(strcmpi({SM.Source}, SurveySources{N})==1)),...
        datestr(SM(LastSurvey).Datenum))
    end
end


