

if exist('/project/group','dir') % iwa
 mpath='/project/group/MOPS/';
else
% set cpg mop folder name
 mpath='/Users/William/Desktop/MOPS/'; % local mac
end

% mpath='/Volumes/group/MOPS/';
% if ~exist(mpath,'dir')
%     fprintf('Using local CGP MOPS folder.\n')
%     mpath='/Users/William/Desktop/MOPS/';
% end

% add matlab paths
eval(['addpath ' mpath]);
eval(['addpath ' mpath 'toolbox']);