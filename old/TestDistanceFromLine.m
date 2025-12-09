% example code to find data points within an absolute distance
% tolerance dtol from a line, where all the x,y data falls between 0 and
% 1 on both axes.

% set the distance tolerance bettwen 0 and 1

dtol=0.1;

% define line y=ax+b
a=2;b=-.2;

% plot the line between x=0 and 1;
x=0:0.01:1;
y=a.*x+b;

figure;
plot(x,y,'k-')
hold on


% generate 100 random x,y points between 0,1
xr=rand(100,1);
yr=rand(100,1);
plot(xr,yr,'k.');
hold on;

% get the distance of each x,y point from the line 
d=DistanceFromLine(xr,yr,a,b);

idx=find(abs(d) <= dtol); % find the points within dtol of line
plot(xr(idx),yr(idx),'ro');
set(gca,'xlim',[0 1],'ylim',[0 1]);


function d=DistanceFromLine(x,y,a,b)

% perpendicular distance of the point x,y from the line defined
%   by y = ax + b

d=abs(a.*x - y + b)./sqrt(a.^2+1);

end