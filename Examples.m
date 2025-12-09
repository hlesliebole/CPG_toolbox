
% find the CA MopNumber (1= Mx border ; 11594= Or border)
%   of the tranditional CountyMopNumber number name.
load MopTableUTM.mat
MopNumber=find(strcmp('L0904',table2cell(Mop(:,1)))); % Malibu

MopFile='Mop00589Data.mat';
load(MopFile,'MopSurvey')

M=load('Mop00582Data.mat','MopSurvey');
M=[M load('Mop00583Data.mat','MopSurvey')];

M(1)=load('Mop00582Data.mat','MopSurvey');
M(2)=load('Mop00583Data.mat','MopSurvey');


fprintf('\n %s contains %s surveys.\n\n',MopFile,num2str(size(MopSurvey,2)))

% isolate individual data sets by their data set number
DataSets=unique(vertcat(MopSurvey.DataSetNum));
fprintf('Database Search Data Set Numbers: %s \n\n',strjoin(cellstr(num2str((DataSets')))))

SurveySources=unique({MopSurvey.source},'stable'); % unique survey sources
fprintf('Survey Sources: %s \n\n',strjoin(SurveySources,', '))

% Example to get indexes of all Gps surveys
GpsIndexes=find(strcmpi({MopSurvey.source}, 'gps')==1);
% Example to get indexes of Gps Jumbo surveys
JumboIndexes=find(contains({MopSurvey.file}, 'jumbo','IgnoreCase',true)==1); 

% Make a table of Survey Source Info
fprintf('%10s %15s %15s \n','Source','# of Surveys','Last Survey')
fprintf('%10s %15s %15s \n','------','------------','-----------')
for N=1:length(SurveySources)
    LastSurvey=find(strcmpi({MopSurvey.source}, SurveySources{N})==1,1,'last');
    if strcmpi(SurveySources{N},'gps')
        LastJumbo=find( contains({MopSurvey.file}, 'jumbo','IgnoreCase',true) == 1,1,'last');
        fprintf('%10s %10i %20s  Jumbos: %10i %20s\n',SurveySources{N},...
        length(find(strcmpi({MopSurvey.source}, SurveySources{N})==1)),...
        datestr(MopSurvey(LastSurvey).datenum),...
        length(find(contains({MopSurvey.file}, 'jumbo',...
          'IgnoreCase',true)==1)),...
          datestr(MopSurvey(LastJumbo).datenum))
    else
        fprintf('%10s %10i %20s\n',SurveySources{N},...
        length(find(strcmpi({MopSurvey.source}, SurveySources{N})==1)),...
        datestr(MopSurvey(LastSurvey).datenum))
    end
end

% Example 2d and 3d plots of data used in gridding

% Example gridding results.
%
% Lidar data in tif files
%  has already been gridded to 1m resolution, but regrid here
%  to fill any gaps in the qc'd data

