close all
mn=month(datetime([SA.Datenum],'convertfrom','datenum')); % month of surveys

nn=28;
mm=find(~isnan(SA(nn).QC));
% XYZ=horzcat(SA(nn).X(mm),...
%     SA(nn).Y(mm),SA(nn).QC(mm));
jn=find(mn == 1);
figure;scatter3(vertcat(SA(jn).X),...
    vertcat(SA(jn).Y),vertcat(SA(jn).Z),'k.');
figure
XYZ=horzcat(vertcat(SA(jn).X),...
    vertcat(SA(jn).Y),vertcat(SA(jn).Z)*0);

%% Step 1. Center the data at zero:
xyz0=mean(XYZ);
A=bsxfun(@minus,XYZ,xyz0); %center the data
scatter3(A(:,1),A(:,2),A(:,3),'k.');
%% Step 2. Find the direction of most variance using SVD and rotate the data to make that the x axis.
[U,S,V]=svd(A,0);
A_rot = A*V; %V(:,1) is the direction of most variance
A_rot(:,3)=vertcat(SA(jn).Z);
hold on, scatter3(A_rot(:,1),A_rot(:,2),vertcat(SA(jn).Z),'k.');

% ibad=isoutlier(A_rot(:,3),'median');%'quartiles');
%   scatter3(A_rot(ibad,1),A_rot(ibad,2),A_rot(ibad,3),'r.');
  %scatter3(A_rot(ibad,1),A_rot(ibad,2),A_rot(ibad,3),'r.');

xstep=10; % xshore qc step size in meters
[xs,iy]=sort(A_rot(:,1));ys=A_rot(iy,2);dzs=A_rot(iy,3);
figure;plot(xs,dzs,'k.');hold on;
xl=1+floor(xs/xstep);
for l=min(xl):max(xl)
    kdx=find(xl == l);
    ibad=isoutlier(dzs(kdx),'median');
    %plot(xs(kdx(ibad)),dzs(kdx(ibad)),'r.')
end
% ibad=isoutlier(dzs,'movmean',round(numel(dzs)/10));
% plot(xs(ibad),dzs(ibad),'r.')