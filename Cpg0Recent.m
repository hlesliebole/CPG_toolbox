fprintf(' \n----------    Recent Surveys added to CPG database  -----------\n')

load SurveyMasterListWithMops.mat;
{Survey(end:-1:end-10).File}'
fprintf('\n ----------    New Truck Files on reefbreak -----------\n\n')
!ls -lt /volumes/group/lidar/VMZ*/L*2/2024*/Beach_Only/*.tif | head -10