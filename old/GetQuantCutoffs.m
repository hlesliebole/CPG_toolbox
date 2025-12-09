
function [qX,pX,pmin,pmax]=GetQuantCutoffs(X,loops)

% initialize min max outlier probability cutoffs as 0,1
%   (no cutoff points found)
pmin=0;
pmax=1;
Xmin=min(X(:));
Xmax=max(X(:));


% get quantile function
pX=0:0.001:1;
qX=quantile(X(:),pX);

% if the X data contains one or a few particularly extreme outliers
%  the slope break points can be set pretty low/high.  The loop
%  setting allows for multiple loops through the data, where the
%  outliers found in each loop are temporarily reset to the median
%  and the slope breaks are recalculated. This could be refined a bit
%  more to look at how many outliers are identified with a particular
%  pmin,pmax pair and the decision to try another loop could be based
%  on how many points were found.

for loop = 1:loops

%-------------------
% min outlier probabilty cutoff
%-------------------
p=0.0:0.001:.1;
tfmax=0;fac=0;
while tfmax == 0 && fac < 1000
fac=fac+10;
tf=ischange(fac*gradient(quantile(X(:),p)),'linear');
tfmax=max(tf);
end

if tfmax > 0
pmin=p(find(tf == 1,1,'last'));
Xmin=quantile(X(:),pmin);
end

%-------------------
% max outlier probabilty cutoff
%-------------------

p=.90:0.001:1;
tfmax=0;fac=0;
while tfmax == 0 && fac < 1000
fac=fac+10;
tf=ischange(fac*gradient(quantile(X(:),p)),'linear');
tfmax=max(tf);
end

if tfmax > 0
pmax=p(find(tf == 1,1,'first'));
Xmax=quantile(X(:),pmax);
end

% set found outliers to the median X within the function 
% and loop through the data more than once to try and get more
% outliers by resetting pmin,pmax.
X(X < Xmin)=nanmedian(X(:));
X(X > Xmax)=nanmedian(X(:));

end

end