% Example code to look at selected MOP CPG process survey files
%  for a MOP number

% Uses the m-script PlotCMstruct.m in MOPS/toolbox

% set path to CPG MOP files from your machine
addpath /Volumes/group/MOPS  % folder with MOP mat files
addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts

% load the "Survey Morphology" *SG.matfile for 582 containing grids
%  for all the surveys
matfile='M00582SG.mat';
load(matfile,'SG')
SG

% list the SG survey sources
fprintf('\n%s contains the following SG.Source survey types:\n',matfile)
unique({SG.Source})'

% list any dates with wheelie surveys

% wheelie survey SG struct arrray index numbers
idx=find(strcmp({SG.Source},'iG8wheel'));

if(numel(idx) > 0)
    fprintf('\n%s contains %i iG8wheel surveys:\n',matfile,numel(idx)')
    fprintf('\nSG Index    Survey Date\n-----------------------\n')
    for n=idx
        fprintf('%8i  %s\n',n,datestr(SG(n).Datenum))
    end
end
    
% plot the most recent survey (last ig8wheel index number)
figure;
PlotCMstruct(SG,idx(end));
title([{matfile};{datestr(SG(idx(end)).Datenum)}]);



