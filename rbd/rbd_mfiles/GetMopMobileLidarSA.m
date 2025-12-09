% Get the SA and QC struct array data for a Mop and identify
%  the mobile lidar surveys in the struct arrays in the index
%   vector "ldx"

% Keep the date of the CurrentSurveyNumber in the current SA array
%  to match with the new SA array.

if exist('SA','var')
  CurrentDatenum=SA(ldx(CurrentSurveyNumber)).Datenum;
else    
  CurrentSurveyNumber=1;
end

% load new SA and QC data
SAmatfile=['M' num2str(CurrentMopNumber,'%5.5i') 'SA.mat'];
fprintf('Loading SA struct array from %s\n',SAmatfile)
load(SAmatfile);
% extract the matlab filepath to the SA.mat file for use when looking for
% and/or saving the QC.mat file
CpgMopFilepath = fileparts(which(SAmatfile));

% full path to the QC mat file
QCmatfile=[CpgMopFilepath filesep 'M' num2str(CurrentMopNumber,'%5.5i') 'QC.mat'];
% get QC struct array if exist, otherwise make one to mirror SA
%  but without any matching x,y,z qc removal points
fprintf('Loading QC struct array from %s\n',QCmatfile)
if exist("QCmatfile","file")
    load(QCmatfile);
else
    fprintf('QC file: %s\n Not found. Making empty QC struct array.\n',QCmatfile)
    QC=SA;
    for n=1:size(QC,2)
        QC(n).X=[]; QC(n).Y=[]; QC(n).Z=[];
    end
end

% identify lidar surveys
ldx=find((strcmp({SA.Source},'AtvMR') | strcmp({SA.Source},'Trk')));
fprintf(' \n    Now Editing MOP  %i \n',CurrentMopNumber);
fprintf(' 1. %i total Mobile LiDAR surveys\n',numel(ldx));
fprintf(' 2. First Survey %s \n',datetime(SA(ldx(1)).Datenum,'convertfrom','datenum'));
fprintf(' 3. Last Survey %s \n',datetime(SA(ldx(end)).Datenum,'convertfrom','datenum'));

% set the CurrentSurveyNumber to as close to the previous Mop current
%  survey datenum as possible.

if CurrentSurveyNumber > 1
    DN=[SA(ldx).Datenum];
    [~,CurrentSurveyNumber]=min(abs(DN-CurrentDatenum));
end

PlotSAQC
edit_menuSA