% make sure group/MOPS/ , group/MOPS/toolbox and group/MOPS/toolbox/profiles %  on reefbreak are in your Matlab path

% load the desired SA struct array
load M00654SA.mat

% now pass to the profile function using a modest YdistTol that should
%  pull in most jumbo data that is near the transect
NumSubTrans=10; % 1 will just return the main transect info; 100 would make
                % statistical profiles with complete use of any 1m res LiDAR
XgapTol=5; % patch any cross-shore profile gaps < 5m wide
YdistTol=25; % can go as large as 50m up- downcoast for a single Mop area

% get nearest point profiles

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);
% (Note: It can take a few minutes to run depending on the number of 
%    NumSubTrans subtransects being considered.)

% take a quick look at all the output Z1Dmedian profiles
figure;
imagesc(X1Dcpg,datenum(Zdatetime),Z1Dmedian,'AlphaData',~isnan(Z1Dmedian));
set(gca,'ydir','normal');datetick('y');
BeachColorbar; % add beach friendly colormap/colorbar from the toolbox
