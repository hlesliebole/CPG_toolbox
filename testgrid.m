load M00474to00495RC.mat
[y,x]=utm2deg(vertcat(RC.Xutm),vertcat(RC.Yutm),...
    repmat('11 S',[length(vertcat(RC.Xutm)) 1]));

% xmin=min(vertcat(RC.Xutm));xmax=max(vertcat(RC.Xutm));
% ymin=min(vertcat(RC.Yutm));ymax=max(vertcat(RC.Yutm));
% xmin=floor(xmin);xmax=ceil(xmax);
% ymin=floor(ymin);ymax=ceil(ymax);
% % x,y grid arrays
% [X,Y]=meshgrid(xmin:xmax,ymin:ymax);
% idx=sub2ind(size(X),vertcat(RC.Yutm)-ymin+1,vertcat(RC.Xutm)-xmin+1);
% Z=X*NaN;Z(idx)=vertcat(RC.Sigma);
% p1=imagesc(X,Y,Z);