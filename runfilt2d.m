
%-----------------------------------------------------------------
function [avgd,stdd]=runfilt2d(wsize,dem,minNR)
%-----------------------------------------------------------------

% returns grids of the average and standard deviation
%  of the input dem array grid point values, for a wsize x wsize
%  indice window centered on each point.

% If a dem grid point was originally a NaN, then the returned
%  average and standard deviations are also returned as NaNs
%  regardless of whether or not there were real values in the
%  surrounding window.  To allow for these points to have real
%  numbers returend instead, comment out the last 2 lines in this
%  function.


% wsize x wsize 2d mean window size
% wsize needs ot be an odd number for convolution so round up 
%   to odd number if necessary
wsize=round(wsize/2)*2+1; 
% make convolution window of ones
W=ones(wsize,wsize);  

% make a mask array dem0 containing only 1's (= real dem numbers) and 0's
% ( = dem NaN values), to figure out the number of real numbers within each
% window-centered grid point, NR.
dem0=dem*0+1; % real numbers to 1's
dem0(isnan(dem0) == 1)=0; % NaNs to 0's
% number of real numbers in each grid point centered window
NR=conv2(dem0,W,'same'); 

% get window avg values
dem(isnan(dem) == 1)=0; % set NaN's to 0
avgd=conv2(dem,W,'same')./NR; % mean window values

% Standard deviation calcs
% Now set avgs at points that are dem 0's to 0's to
% avoid messing up std calcs
avgd0=avgd;avgd0(dem == 0)=0; 
asq=conv2(dem.^2,W,'same'); % square of dem values
bsq=conv2(avgd0.^2,W,'same'); % square of mean values
m2ab=conv2(2.*dem.*avgd0,W,'same');  % cross term in stdev equation
stdd=sqrt((asq-m2ab+bsq)./NR); % window standard deviations

% Reset grid points that were originally NaNs back to NaN's based
%  on zero points in dem0.  Comment these lines out if you want
%  to return any potential real number results at these points
%  (ie there were real numbers in the window centered on the NaN point)
avgd(dem0 == 0)=NaN;
stdd(dem0 == 0)=NaN;

% only keep window results with > 95% coverage (drops edges
%  mostly)
% avgd(NR < 0.95*(wsize.^2))=NaN;
% stdd(NR < 0.95*(wsize.^2))=NaN;
avgd(NR < minNR)=NaN;
stdd(NR < minNR)=NaN;

end
