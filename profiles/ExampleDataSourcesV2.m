
MopFile=['M00654SA.mat'];
load(MopFile)
fprintf('\n %s contains %s surveys.\n\n',MopFile,num2str(size(SA,2)))

SurveySources=unique({SA.Source},'stable'); % unique survey sources
fprintf('Survey Sources: %s \n\n',strjoin(SurveySources,', '))

% Example to get indexes of all Gps surveys
GpsIndexes=find(strcmpi({SA.Source}, 'gps')==1);
% Example to get indexes of Gps Jumbo surveys
JumboIndexes=find(contains({SA.File}, 'jumbo','IgnoreCase',true)==1); 

% Make a table of Survey Source Info
fprintf('%10s %15s %15s %15s\n','Source','# of Surveys','First Survey','Last Survey')
fprintf('%10s %15s %15s %15s\n','------','------------','------------','-----------')

for N=1:length(SurveySources)

    FirstSurvey=find(strcmpi({SA.Source}, SurveySources{N})==1,1,'first');
    LastSurvey=find(strcmpi({SA.Source}, SurveySources{N})==1,1,'last');
    FirstJumbo=find( contains({SA.File}, 'jumbo','IgnoreCase',true) == 1,1,'first');
    LastJumbo=find( contains({SA.File}, 'jumbo','IgnoreCase',true) == 1,1,'last');
    
    if strcmpi(SurveySources{N},'gps')
    
        fprintf('%10s %10i % 20s %15s\n',[SurveySources{N} ':All'],...
          length(find(strcmpi({SA.Source}, SurveySources{N})==1)),...
          datestr(SA(FirstSurvey).Datenum),datestr(SA(LastSurvey).Datenum))
   
        fprintf('%10s %10i %20s %15s\n',[SurveySources{N} ':Jumbos'],...
          length(find(contains({SA.File}, 'jumbo','IgnoreCase',true)==1)),...
          datestr(SA(FirstJumbo).Datenum),datestr(SA(LastJumbo).Datenum))

    else

        fprintf('%10s %10i %20s %15s\n',SurveySources{N},...
          length(find(strcmpi({SA.Source}, SurveySources{N})==1)),...
          datestr(SA(FirstSurvey).Datenum),datestr(SA(LastSurvey).Datenum))
    
    end

end


