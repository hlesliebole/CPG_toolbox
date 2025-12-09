% example code to view all the shoreline struct array data

% load shore mat file
load TorreyPinesShoreline.mat
MopStart=min(SS(1).MopNumbers);
MopEnd=max(SS(1).MopNumbers);

% load Matlab table of Mop transects info
load('MopTableUTM.mat','Mop')

%---------------------------------------------------------------
% plot all the shoreline locations as different colored lines

%---------------  MSL shorelines -----------------------------------------
% FIGURE 1
figure('position',[14  146  1060  651],'name','MSL Shorelines'); 
subplot(1,2,1);hold on;
for i=1:size(SS,2);plot(SS(i).MSLlongitudes,SS(i).MSLlatitudes,'-');end
for m=MopStart:MopEnd;pt=plot([Mop.BackLon(m) Mop.OffLon(m)],[Mop.BackLat(m) Mop.OffLat(m)],'k-');end 
p1=plot(GS.MSLlongitudes,GS.MSLlatitudes,'m-','linewidth',2);
legend([p1 pt],'Global Mean MSL Shoreline','Mop Transects 520-596','location','southwest');
xlabel('Longitude');ylabel('Latitude');title('All Surveys MSL Lat-Lon Shorelines'); 
subplot(1,2,2);hold on;
for i=1:size(SS,2);plot(SS(i).MSLeastings,SS(i).MSLnorthings,'-');end
for m=MopStart:MopEnd;pt=plot([Mop.BackXutm(m) Mop.OffXutm(m)],[Mop.BackYutm(m) Mop.OffYutm(m)],'k-');end
p1=plot(GS.MSLeastings,GS.MSLnorthings,'m-','linewidth',2);
legend([p1 pt],'Global Mean MSL Shoreline','Mop Transects 520-596','location','southwest');
xlabel('Eastings');ylabel('Northings');title('All Surveys MSL UTM Shorelines'); 

%------  MHW shorelines ----------------------------------------
% FIGURE 2
figure('position',[44  126  1060  651],'name','MHW Shorelines');
subplot(1,2,1);hold on;
for i=1:size(SS,2);plot(SS(i).MHWlongitudes,SS(i).MHWlatitudes,'-');end
for m=MopStart:MopEnd;pt=plot([Mop.BackLon(m) Mop.OffLon(m)],[Mop.BackLat(m) Mop.OffLat(m)],'k-');end 
p1=plot(GS.MHWlongitudes,GS.MHWlatitudes,'m-','linewidth',2);
legend([p1 pt],'Global Mean MHW Shoreline','Mop Transects 520-596','location','southwest');
xlabel('Longitude');ylabel('Latitude');title('All Surveys MHW Lat-Lon Shorelines'); 
subplot(1,2,2);hold on;
for i=1:size(SS,2);plot(SS(i).MHWeastings,SS(i).MHWnorthings,'-');end
for m=MopStart:MopEnd;pt=plot([Mop.BackXutm(m) Mop.OffXutm(m)],[Mop.BackYutm(m) Mop.OffYutm(m)],'k-');end
p1=plot(GS.MHWeastings,GS.MHWnorthings,'m-','linewidth',2);
legend([p1 pt],'Global Mean MHW Shoreline','Mop Transects 520-596','location','southwest');
xlabel('Eastings');ylabel('Northings');title('All Surveys MHW UTM Shorelines'); 

% plot beach widths and shoreline location anomalies as a function
%  of the Mop number
%-----  MSL beach widths and shoreline location anomalies -----------------------------------------
% FIGURE 3
figure('position',[74  106  1060  651],'name','MSL Beach Widths and Shoreline Location Anomalies'); 
subplot(2,1,1);hold on;
for i=1:size(SS,2);plot(SS(i).MopNumbers,SS(i).MSLbeachWidths,'-');end
p1=plot(GS.MopNumbers,GS.MSLbeachWidths,'m-','linewidth',2);set(gca,'xlim',[520 596]);
legend(p1,'Global Mean MSL Beach Width','location','northwest');grid on;
xlabel('Mop Number');ylabel('Beach Width (m)');title('All Surveys MSL Beach Width'); 
subplot(2,1,2);hold on;
for i=1:size(SS,2);plot(SS(i).MopNumbers,SS(i).MSLanomalies,'-');end
set(gca,'xlim',[520 596]);grid on;
xlabel('Mop Number');ylabel('Location Anomaly (m)');
title('All Surveys MSL Shoreline Location Anomalies relative to Global Mean Shoreline'); 
%---------------------------------------------------------------
% FIGURE 4
figure('position',[104  86  1060  651],'name','MHW Beach Widths and Shoreline Location Anomalies'); 
subplot(2,1,1);hold on;
for i=1:size(SS,2);plot(SS(i).MopNumbers,SS(i).MHWbeachWidths,'-');end
p1=plot(GS.MopNumbers,GS.MHWbeachWidths,'m-','linewidth',2);set(gca,'xlim',[520 596]);
legend(p1,'Global Mean MHW Beach Width','location','northwest');grid on;
xlabel('Mop Number');ylabel('Beach Width (m)');title('All Surveys MHW Beach Width'); 
subplot(2,1,2);hold on;
for i=1:size(SS,2);plot(SS(i).MopNumbers,SS(i).MHWanomalies,'-');end
set(gca,'xlim',[520 596]);grid on;
xlabel('Mop Number');ylabel('Location Anomaly (m)');
title('All Surveys MHW Shoreline Location Anomalies relative to Global Mean Shoreline'); 

%---------------------------------------------------------------
% 3D color scatter plot of all the MSL-MHW shoreface slopes
% FIGURE 5
figure('position',[134  66  1060  651],'name','Mop  MSL-MHW shoreface slopes');
x=repmat([SS.Datenum],77,1);x=x(:)'; % make survey date for every points
y=[SS.MopNumbers];
z=(1.344-0.774)./([SS.MSLbeachWidths]-[SS.MHWbeachWidths]);
amax=0.15; % set max slope for coloring
zscaled = 1+64*(z)/(amax);
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled) & z > 0); % non NaN points and positive slopes
                                     
cm =jet(64);

scp=scatter3(x(idx), y(idx), z(idx), 15, cm(ceil(zscaled(idx)),:), 'filled');
scp.MarkerFaceAlpha = .9;
scp.MarkerEdgeAlpha = .9;
view(-10,80)
colormap(jet(64))
cb=colorbar;cb.Label.String='Shoreface Slope';
set(gca,'clim',[0 amax]);
set(gca,'xlim',[min(x) max(x)]);
set(gca,'zlim',[min(z) max(z)]);
set(gca,'color',[.7 .7 .7]);
datetick;
xlabel('Survey Date');ylabel('Mop Number');zlabel('Anomaly (m)');
title('All Surveys MSL-MHW Shoreface Slopes'); 



%---------------------------------------------------------------
% 3D color scatter plot of all the MSL shoreline anomalies
% FIGURE 6
figure('position',[164  46  1060  651],'name','Mop Shoreline Location Anomalies');
x=repmat([SS.Datenum],77,1);x=x(:)'; % make survey date for every points
y=[SS.MopNumbers];
z=[SS.MSLmonthanomalies];
amax=ceil(max(abs(z)));
amax=56;
zscaled = 1+64*(z+amax)/(2*amax);
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points
                                     
cm =flipud(polarmap(64));

scp=scatter3(x(idx), y(idx), z(idx), 15, cm(ceil(zscaled(idx)),:), 'filled');
scp.MarkerFaceAlpha = .9;
scp.MarkerEdgeAlpha = .9;
view(-10,80)
colormap(flipud(polarmap(64)))
cb=colorbar;cb.Label.String='Shoreline Location Anomaly (m)';
set(gca,'clim',[-amax amax]);
set(gca,'xlim',[min(x) max(x)]);
set(gca,'zlim',[min(z) max(z)]);
set(gca,'color',[.7 .7 .7]);
datetick;
xlabel('Survey Date');ylabel('Mop Number');zlabel('Anomaly (m)');
title('All Surveys MSL Shoreline Location Anomalies relative to Global Mean Shoreline'); 



% 
% MSL2d=reshape([SS.MSLanomalies],[length(SS(1).MSLbeachWidths),size(SS,2)]);
% pcolor([SS.Datenum],[SS(1).MopNumbers],MSL2d);datetick;colormap(jet);colorbar;
% xlabel('Date');ylabel('Mop Number');shading flat;
% 
% figure
% for i=1:size(SS,2)
%     %plot(SS(i).MHWeastings,SS(i).MHWnorthings,'-')
%     plot(SS(i).MHWlongitudes,SS(i).MHWlatitudes,'-')
%     hold on;
% end
% for m=MopStart:MopEnd
%     %plot([Mop.BackXutm(m) Mop.OffXutm(m)],[Mop.BackYutm(m) Mop.OffYutm(m)],'k-');
%     plot([Mop.BackLon(m) Mop.OffLon(m)],[Mop.BackLat(m) Mop.OffLat(m)],'k-'); 
% end
% plot(GS.MHWlongitudes,GS.MHWlatitudes,'m-','linewidth',2)
% %xlabel('UTM Easting');ylabel('UTM Northing');
% xlabel('Longitude');ylabel('Latitude');
%   

function cmap = polarmap(varargin)
%POLARMAP Polarized color map
%	POLARMAP applies a "polarized" blue-white-red colormap to current figure,
%	and adjusts the color axis limits to be centered to zero.
%
%	POLARMAP(M) fixes the number of colors to M (default is 64).
%
%	POLARMAP(MAP) applies linear shading to white to the center of colormap
%	MAP which can be any of existing colormaps (an Mx3 matrix of RGB).
%
%	POLARMAP(MAP,C) uses exponent C to modify the shading contrast. Default 
%	is C = 1 for linear shading. Use C = 2 to strengthen the shading, or 
%	C = 0.5 to attenuate it.
%
%	C=POLARMAP(...) returns an M-by-3 matrix containing the colormap, that 
%	can be used with COLORMAP function like other colormaps.
%
%	Examples:
%		pcolor(peaks), shading interp
%		polarmap, colorbar
%
%	then try the following
%		polarmap(jet,0.5)
%
%	Note the polar shading has no real interest with colormaps that include
%	white color as one of the extremes (like GRAY, BONE, HOT, ...).
%
%	See also JET, HSV, COPPER, SPRING, SUMMER, WINTER, COOL, COLORMAP, RGBPLOT.
%
%	Author: Francois Beauducel, IPGP
%	Created: 2011-10-26
%	Updated: 2012-06-12

%	Copyright (c) 2012, François Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without 
%	modification, are permitted provided that the following conditions are 
%	met:
%
%	   * Redistributions of source code must retain the above copyright 
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright 
%	     notice, this list of conditions and the following disclaimer in 
%	     the documentation and/or other materials provided with the distribution
%	                           
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%	POSSIBILITY OF SUCH DAMAGE.

% default parameters
m = 64;	% number of colors
c = 1;	% exponent of shading factor (1 = linear)

if nargin > 0
	if ~isnumeric(varargin{1}) | (size(varargin{1},2) ~= 3 & ~isscalar(varargin{1}))
		error('First argument must be numeric: scalar M or Mx3 color matrix');
	end
	if isscalar(varargin{1})
		m = varargin{1};
	end
end
if nargin > 0 & size(varargin{1},2) == 3
		map = varargin{1};
		m = size(map,1);
else
	map = bluered(m);
end

if nargin > 1 & isscalar(varargin{2})
	c = varargin{2};
end

% linear shading from min/max (colormap value) to center (white)
r = repmat(abs(linspace(1,-1,m)).^c,[3,1])';
map = map.*r + 1 - r;

if nargout > 0
	cmap = map;
else
	colormap(map)
	caxis([-1,1]*max(abs(caxis)))
	% Note: this fixes color axis to manual mode...
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map = bluered(m)

if mod(m,2)
	z = [0,0,0];
	m2 = floor(m/2);
else
	z = zeros([0,3]);
	m2 = m/2;
end
map = [repmat([0,0,1],[m2,1]);z;repmat([1,0,0],[m2,1])];
end