% example script to calculate a mop xshore profile in
% Mop transect x,z coordinates from the 1m averaged 
% survey points in the SA cpg mop file

addpath /volumes/group/MOPS/toolbox
addpath /volumes/group/MOPS

% load SA struct array in mop 582 SA mat file
load M00582SA.mat

% find the index of the jumbo on 12 6 2019
idx=find([SA.Datenum] == datenum(2019,12,6) & ...
    contains({SA.File}, 'jumbo','IgnoreCase',true) );
    
%  convert SA 1m avg survey point x,y utm values
%  to cross-shore transect X distances in whole meters.
%  This generic function for any utm x,y points returns
%  arrays of the closest mop numbers (Nmop) and their 
%  mop xshore X values for each point.  In this application,
%  using SA struct array data, all the x,y data is for the same 
%  mop so all the returned Nmop values are the same = 582.
[Nmop,X]=UTM2MopxshoreX(SA(idx).X,SA(idx).Y); 

% now average the survey X,z values into unique 1m xshore X distance bins
[ux, ~, xidx] = unique(X); % ux are unique 1m resolution offshore distances
zavg=accumarray(xidx(:), SA(idx).Z',[], @mean); % mean z for each ux

% now interpolate in ux to make sure to have a continuous profile from 
%  the min to max ux values.  
X1D=min(ux):max(ux);
Z1D=interp1(ux,zavg,X1D);

% make plot
figure;plot(X1D,Z1D,'k-');
set(gca,'xdir','reverse');xlabel('Xshore Distance (m)');ylabel('z (m.navd88)');
title({'xshore-interpolated average profile',...
    'from 1m spatial average survey data points in SA struct array',...
    datestr(SA(idx).Datenum,'mm/dd/YYYY')})


