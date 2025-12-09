function [X,Y,Z]=SG2grid(SG,SurvNum)

MopNumber=SG(1).Mopnum;

% Set generic grid limits based on all survey grid data for this mop.
% This makes it easy to make a difference grid from mutiple calls of
% this function for different surveys.

xmin=min(vertcat(SG.X));xmax=max(vertcat(SG.X));
ymin=min(vertcat(SG.Y));ymax=max(vertcat(SG.Y));

% xmin=min([Mop.BackXutm(MopNumber-1:MopNumber+1)' Mop.OffXutm(MopNumber-1:MopNumber+1)']);
% xmax=max([Mop.BackXutm(MopNumber-1:MopNumber+1)' Mop.OffXutm(MopNumber-1:MopNumber+1)']);
% ymin=min([Mop.BackYutm(MopNumber-1:MopNumber+1)' Mop.OffYutm(MopNumber-1:MopNumber+1)']);
% ymax=max([Mop.BackYutm(MopNumber-1:MopNumber+1)' Mop.OffYutm(MopNumber-1:MopNumber+1)']);
% xmin=min([xmin [SG(SurvNum).X]']);xmax=max([xmax [SG(SurvNum).X]']);
% ymin=min([ymin [SG(SurvNum).Y]']);ymax=max([ymax [SG(SurvNum).Y]']);

xmin=floor(xmin);xmax=ceil(xmax);
ymin=floor(ymin);ymax=ceil(ymax);

x=SG(SurvNum).X;
y=SG(SurvNum).Y;

% x,y grid arrays
[X,Y]=meshgrid(xmin:xmax,ymin:ymax);
idx=sub2ind(size(X),y-ymin+1,x-xmin+1);

% z grid array
Z=X*NaN; % initialize as NaNs
Z(idx)=SG(SurvNum).Z; % assign valid z grid data points

% figure
% p1=surf(X,Y,Z,'linestyle','none');    
% grid on;
% set(gca,'color',[.7 .7 .7],'fontsize',14);
% y_labels = get(gca, 'YTick');
% set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
% ylabel('northings (m)');
% x_labels = get(gca, 'XTick');
% set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
% xlabel('eastings (m)');

end
