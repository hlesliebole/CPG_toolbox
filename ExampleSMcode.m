% Example code to look at selected MOP CPG processed survey profiles
%  for a single MOP number

% set path to CPG MOP files from your machine (eg. for a mac)
addpath /Volumes/group/MOPS  % folder with MOP mat files
addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts

% load the "Survey Morphology" *SM.matfile for TP Mop 582 containing profiles
%  for all the surveys
matfile='M00644SM.mat';
load(matfile,'SM')
SM

% list the SM survey sources
fprintf('\n%s contains the following SM.Source survey types:\n',matfile)
unique({SM.Source})'

% list any dates with wheelie surveys

% truck lidar survey SM struct arrray index numbers
idx=find(strcmp({SM.Source},'Trk'));

if(numel(idx) > 0)
    fprintf('\n%s contains %i Trk surveys. Plotting the most recent one.\n',matfile,numel(idx)')
    % fprintf('\nSM Index    Survey Date\n-----------------------\n')
    % for n=idx
    %     fprintf('%8i  %s\n',n,datestr(SM(n).Datenum))
    % end
end
    
% plot the most recent truck survey (last Trk index number)
figure;
PlotCMstruct(SM,idx(end));
title([{matfile};{datestr(SM(idx(end)).Datenum)}]);



