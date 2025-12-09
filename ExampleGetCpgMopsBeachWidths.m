% Example code to calculate a time series of Mop beach widths
%  relative to the Mop back beach point for historic survey data 
%  in the CPG MOP database.

%% Settings
% path settings (change these if you are using a windows machine)
addpath /volumes/group/Mops   % location of Mop survey database
addpath /volumes/group/Mops/toolbox  % location of supporting m-scripts 

% mop number (sequential numeric number from MX border convention. For
% San Diego county, this is the same as the D**** Mop number without the D.)
MopNumber=580; % example for Torrey Pines parking lot point (lots of survey data)

%% Load the processed profile info in the Survey Morphology "SM"
%  mat file for the mop number
load(['M00' num2str(MopNumber,'%2.2i') 'SM.mat'])

%% Loop through all the survey and find the profile intersections with
%  MSL and MHW, if they exist. 

Sdate=NaN(size(SM,2),1); % preallocate survey date vector
Xmsl=NaN(size(SM,2),1); % preallocate MSL beach width as no data NaN
Xmhw=NaN(size(SM,2),1); % preallocate MHW beach width as no data NaN
for n=1:size(SM,2)   
        
        Sdate(n)=SM(n).Datenum; % survey datenumber
       
        % Use Mops/toolbox/interesections.m to find widths
        %  In this example, I'm using the profile directly along the Mop
        %  line (.Z1Dtransect) interpolated from the gridded survey data.
        %  You can also use .Z1Dmean (mean profile for Mop area) and
        %  .Z1Dmedian (median profile for the mop area).
        xMSL=intersections([-10 150],[.774 .774],SM(n).X1D,SM(n).Z1Dtransect);
        xMHW=intersections([-10 150],[1.344 1.344],SM(n).X1D,SM(n).Z1Dtransect);
        
        % Sometimes there can be no intersections, so leave value as NaN,
        % or mutiple intersections (use the one closest to the back beach
        % point here, but you might prefer to do something else) - or
        % print a warning when more than 1 intersection is found etc.
        
        if numel(xMSL) > 0 % screen for any intersections
          Xmsl(n)=min(xMSL); % takes smallest width if more than one
        end
        
        if numel(xMHW) > 0 % screen for any intersections
          Xmhw(n)=min(xMHW); % takes smallest width if more than one
        end
end

%  plot the results
figure;
plot(Sdate,Xmsl,'b-','DisplayName','MSL');
hold on;
plot(Sdate,Xmhw,'s','DisplayName','MHW');
legend
datetick
grid on
xlabel('Date');ylabel('Beach Width (m) from Mop Back Beach Point');
title(['Mop:' num2str(MopNumber)]);


