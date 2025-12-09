function [X0BeachOnlyTrk,X0BeachOnlyAtvMR,nTrkSurveys,nAtvSurveys]=GetBeachOnlyX0(SA)

% Finds the mobile LiDAR back beach location along a mop line, X0,
% in integer meters along the transect relative to the officeal Mop back beach
% point location.  Returns separate values based on the truck LiDAR
% surveys (preferred choice) and the AtvMR mini-ranger LiDAR surveys
% as these are processed and QC'd somewhat separately.  X0 is the median
% last transect line point with valid LiDAR data for all the Trk or AtvMR
% surveys.

% input: 
%        SA struct array for a Mop

% output:
%         X0BeachOnlyTrk = Trk LiDAR data Mop transect back beach point 
%                           distance from Mop back beach point (m)
%
%         X0BeachOnlyAtvMR = AtvMR LiDAR data Mop transect back beach point 
%                           distance from Mop back beach point (m)

%% ------------------------------------------------

nTrkSurveys=0;nAtvSurveys=0;
%% reduce SA to just Trk and AtvMR surveys
idx=find(strcmp({SA.Source},'Trk') | strcmp({SA.Source},'AtvMR'));
SA=SA(idx);

%% check if thes mop has any mobile LiDAR data
if numel(idx) == 0
    fprintf('No Trk or AtvMR surveys. Setting both X0 = 0m\n')
    X0BeachOnlyTrk=0;
    X0BeachOnlyAtvMR=0;
else
%% find median survey location of back beach point

%% get nearest point profiles with tight mobile LiDAR tolerances
    Ytol=5; % 5m alongcoast nearest point tolerance
    Xtol=3; % 3m cross shore gap interpolation tolerance
    [X1D,Z1D]=GetAllNearestPointsProfiles(SA,Ytol,Xtol);
 
%% find the truck back beach boundary along the mop tranect
    trk=find(strcmp({SA.Source},'Trk'));
    if ~isempty(trk)
    for n=1:numel(trk)
        xback(n)=X1D(find(~isnan(Z1D(trk(n),:)), 1 ));
    end
    else
        xback=0;
        fprintf('No Trk surveys. Setting X0BeachOnlyTrk = 0m\n')
    end
    X0BeachOnlyTrk=median(xback);
    nTrkSurveys=numel(find(~isnan(xback)));
    fprintf('Mop %5i, Median Truck Back Beach X0= %5im, Nsurveys= %4i, min-max X0 = %5im to %5im\n',...
        SA(1).Mopnum,X0BeachOnlyTrk,nTrkSurveys,min(xback),max(xback))

%% find the AtvMR back beach boundary along the mop tranect
    xback=[];
    atv=find(strcmp({SA.Source},'AtvMR'));
    if ~isempty(atv)
    for n=1:numel(atv)
        xback(n)=X1D(find(~isnan(Z1D(atv(n),:)), 1 ));
    end
    else
        xback=0;
        fprintf('No AtvMR surveys. Setting X0BeachOnlyAtvMR = 0m\n')
    end
    X0BeachOnlyAtvMR=median(xback);
    nAtvSurveys=numel(find(~isnan(xback)));
    fprintf('Mop %5i, Median AtvMR Back Beach X0= %5im, Nsurveys= %4i, min-max X0 = %5im to %5im\n',...
        SA(1).Mopnum,X0BeachOnlyAtvMR,nAtvSurveys,min(xback),max(xback))

    
end  % end if-then mobile lidar data exists

end