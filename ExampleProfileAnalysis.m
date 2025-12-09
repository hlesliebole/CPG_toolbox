clearvars
Ytol=50; % 50m alongcoast nearest point tolerance
Xtol=10; % 5m cross shore gap interpolation tolerance
%for Mop=1:746

load M00667SA.mat;
nsurv=size(SA,2);
fprintf('%i surveys.\n',nsurv)
[X1D,Z1D]=GetNearestPointsProfiles(SA,Ytol,Xtol);
X2D=repmat(X1D,size(Z1D,1),1);
trk=find(strcmp({SA.Source},'Trk'));
%Z1D=Z1D(trk,:);
for n=1:numel(trk)
    xback(n)=X1D(find(~isnan(Z1D(trk(n),:)), 1 ));
end
X0BeachOnly=median(xback);
fprintf('Truck Back Beach x= %6.1f \n',X0BeachOnly)
% reduce to common back beach x location based on truck beach-only data
idx=find(X1D >= X0BeachOnly);
Z1Dt=Z1D(:,idx);
X1Dt=X1D(idx);

% figure
% surf(X1D,1:size(Z1D(trk,:),1),Z1D(trk,:));shading flat;BeachColorbar;view(2)
% 
% figure
% surf(X1D,1:size(Z1D,1),Z1D);shading flat;BeachColorbar;view(2)
% 
% figure;contour(X1D,1:size(Z1D,1),Z1D,[1.344 1.344],'r-');
% hold on;
% contour(X1D,1:size(Z1D,1),Z1D,[0.774 0.774],'b-');

%zc=-0.058; %MLLW
%zc=0.218; %MLW
zc=0.774;%MSL
zc=1.344; %MHW
%zc=1.566; %MHHW
%zc=2.138; %HAT

% remove any elevations seaward of the min profile elevation
%  (lidar swash filter)



Zf=Z1Dt;
[zmin,imin]=min(Zf');
for n=1:size(Zf,1)
    if imin(n) > 0
        Zf(n,imin(n)+1:end)=NaN;
    end
end

X2Df=repmat(X1Dt,size(Zf,1),1);

% 
% 
% Z=Zf;
% Z(Z > zc)=NaN;
% [zmax,imax]=max(Z');
% imax(isnan(zmax))=NaN;
% imax(zc-zmax > 0.1)=NaN;
% zmax(zc-zmax > 0.1)=NaN;
% figure;
% subplot(2,1,1);plot(zmax);hold on;subplot(2,1,2);plot(imax);hold on;plot(imax,'o')
% 
% Z=Zf;
% Z(Z < zc)=NaN;
% [zmin,imin]=min(Z');
% imin(isnan(zmin))=NaN;
% imin(zmin-zc > 0.1)=NaN;
% zmin(zmin-zc > 0.1)=NaN;
% %figure;
% subplot(2,1,1);plot(zmin);subplot(2,1,2);plot(imin);hold on;plot(imin,'^')
% imean=(nanmean(vertcat(imin,imax)));plot(imean,'+')
% fprintf('%i surveys have valid %7.3f m contour location data.\n',numel(find(~isnan(imean))),zc)
dztol=0.1;
Z=Zf;
%Z(Z < zc)=NaN;
[zmin,imin]=min(abs(Z-zc)');
imin(isnan(zmin))=NaN;
imin(zmin > dztol)=NaN;
%zmin(zmin-zc > 0.1)=NaN;
figure
subplot(2,1,1);plot(zmin);hold on;plot([0 numel(zmin)],[dztol dztol],'k--');
set(gca,'yscale','log');subplot(2,1,2);plot(imin);hold on;plot(imin,'s')
fprintf('%i surveys have a valid %6.3f m contour location data.\n',numel(find(~isnan(imin))),zc)
% load M00582SM.mat
% Z=reshape([SM.Z1Dtransect],size(SM(1).X1D,2),size(SM,2))';
% 
% TF=isoutlier(Z);
% Zf=Z;Zf(TF == 1)=NaN;
% Zmin=min(Zf);

figure('position',[47          56        1342         739]);
plot(X2D(:),Z1D(:),'k.',X2Df(:),Zf(:),'r.')
hold on;plot(X1D,Z1D(end,:),'g.');
title(SA(end).File)
plot([X0BeachOnly X0BeachOnly],[0 3],'k-');
set(gca,'xdir','reverse')