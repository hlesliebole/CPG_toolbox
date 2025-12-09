% Example code to look at selected MOP CPG process survey files
%  for a MOP number

% Uses the m-script PlotCMstruct.m in MOPS/toolbox

% set path to CPG MOP files from your machine
addpath /Volumes/group/MOPS  % folder with MOP mat files
addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts

% load the "Survey Morphology" *SM.matfile for Mop 38 containing profiles
%  for all the surveys
matfile='M00038SM.mat';
load(matfile,'SM')
SM

% list the SM survey sources
fprintf('\n%s contains the following SM.Source survey types:\n',matfile)
unique({SM.Source})'

% list any dates with wheelie surveys

% wheelie survey SM struct arrray index numbers
idx=find(strcmp({SM.Source},'iG8wheel'));

if(numel(idx) > 0)
    fprintf('\n%s contains %i iG8wheel surveys:\n',matfile,numel(idx)')
    fprintf('\nSM Index    Survey Date\n-----------------------\n')
    for n=idx
        fprintf('%8i  %s\n',n,datestr(SM(n).Datenum))
    end
end
    
% plot the most recent survey (last ig8wheel index number)
figure;
PlotCMstruct(SM,idx(end));
title([{matfile};{datestr(SM(idx(end)).Datenum)}]);



