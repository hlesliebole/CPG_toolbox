
% ExampleQC.m

% Demo of ElevGoodnessQcv1.m using Truck LiDAR 1m spatial average
%  tif x,y,z data that has been sorted into individual Mop areas
%  for LJ Shores.

% pick a Shores Mop number
MopNum=503;

% load the Mop database spatially averaged survey data struct array

load([ 'M'  num2str( MopNum , '%5.5i' )  'SA.mat' ],'SA'); 

% use the 6th survey in the struct array, assign values to x,y,z variables
n=6;
x=SA(n).X;y=SA(n).Y;z=SA(n).Z;

% make before and after QC plots using histogram bin size of 0.01m

figure('position',[30         236        1350         524]);

% show data before qc
subplot(1,2,1);ScatterPlotBeachUTM(x,y,z,'3d');
title({['MOP # ' num2str(MopNum)] , datestr(SA(n).Datenum)})
% qc data
dz=0.01;
g=ElevGoodnessQCv1(x,y,z,dz);

% show QC data when eliminating g=0 points
gthreshold=0;
z2=z;z2(g <= gthreshold)=NaN;subplot(1,2,2);ScatterPlotBeachUTM(x,y,z2,'3d');
title(['QC with Bin Size = ' num2str(dz,'%5.3f')...
    ' m  ; goodness threshold = ' num2str(gthreshold) ])

% make plot of goodness sensitivity to dz
PlotGoodnessVsBinSize
 
 

