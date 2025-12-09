% Updates the SA, SG and SM matfiles for a single survey
% defined the the survey date and the source name
%

SingleDate=datenum(2021,10,22)
SingleSource='Gps'

% source choices:
%     {'Gps'     }
%     {'UTAir'   }
%     {'Trk'     }
%     {'TrkMR'   }
%     {'AtvMR'   }
%     {'DrnMR'   }
%     {'USACE'   }
%     {'CCC'     }
%     {'USGS'    }
%     {'KMair'   }
%     {'iG8wheel'}
%

addpath /volumes/group/mops
addpath /volumes/group/mops/toolbox

cd /volumes/group/mops/toolbox
MakeSurveyMasterList

cd /volumes/group/mops

UpdateSingleSurveySAmatfiles

UpdateSingleSurveySGmatfiles

UpdateSingleSurveySMmatfiles