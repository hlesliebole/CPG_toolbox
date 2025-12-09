function [X0BeachOnly,X1Dt,Zf,Zm,Zm3,Zq,Zs,Zg]=GetMeanMedianProfiles(SM)

%% returns mean profiles for the input SA struct array
%
%  X0BeachOnly = Mop xshore location of truck back beach boundary
%  X1Dt(N) = N profile xshore x values (1m res) relative to X0BeachOnly
%  Zf(M,N) = M survey profiles relative to X0BeachOnly
%  Zm(12,N) = monthly mean profiles (calculated by year first and then across all years)
%  Zq(4,N) = quarterly mean profiles (Jan-Mar,Apr-Jun,Jul-Sep,Oct-Dec)
%  Zs(2,N) = seasonal mean profiles (Jan-Jun,Jul-Dec)
%  Zg(1,N) = global mean profile

%% get nearest point profile
Ytol=50; % 50m alongcoast nearest point tolerance
Xtol=10; % 5m cross shore gap interpolation tolerance

%[X1D,Z1D]=GetNearestPointsProfiles(SA,Ytol,Xtol);
X1D=SM(1).X1D;
Z1D=reshape([SM.Z1Dmedian],size(SM(1).X1D,2),size(SM,2))';

%X2D=repmat(X1D,size(Z1D,1),1);
%% find the truck back beach boundary along the mop tranect
trk=find(strcmp({SM.Source},'Trk'));
if ~isempty(trk)
for n=1:numel(trk)
    xback(n)=X1D(find(~isnan(Z1D(trk(n),:)), 1 ));
end
else
    xback=0;
end
X0BeachOnly=median(xback);
% fprintf('Truck Back Beach x= %6.1f \n',X0BeachOnly)

%% reduce profiles to common back beach x location based on truck beach-only data
idx=find(X1D >= X0BeachOnly);
Z1Dt=Z1D(:,idx);
X1Dt=X1D(idx);

%% remove any elevations seaward of the min profile elevation
%  (lidar swash filter)
Zf=Z1Dt;
[zmin,imin]=min(Zf');
for n=1:size(Zf,1)
    if imin(n) > 0
        Zf(n,imin(n)+1:end)=NaN;
    end
end

% remove any additional outliers
TF=isoutlier(Zf,"mean");
Zf(TF == 1)=NaN;

%X2Df=repmat(X1Dt,size(Zf,1),1);

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
% dztol=0.1;
% Z=Zf;
% %Z(Z < zc)=NaN;
% [zmin,imin]=min(abs(Z-zc)');
% imin(isnan(zmin))=NaN;
% imin(zmin > dztol)=NaN;
% %zmin(zmin-zc > 0.1)=NaN;
% figure
% subplot(2,1,1);plot(zmin);hold on;plot([0 numel(zmin)],[dztol dztol],'k--');
% set(gca,'yscale','log');
% xlabel('Survey Number');ylabel('Dist from contour (m)');
% title(['Min Profile Point Distance from ' num2str(zc) ' m contour']);
% subplot(2,1,2);plot(imin(~isnan(imin)));hold on;plot(X1Dt(imin(~isnan(imin))),'s')
% fprintf('%i surveys have a valid %6.3f m contour location data.\n',numel(find(~isnan(imin))),zc)
% xlabel('Survey Number');ylabel('Dist from back beach (m)');
% title(['Contour Distance from (Truck) back beach (m)' num2str(zc) ' m contour']);
% % load M00582SM.mat
% % Z=reshape([SM.Z1Dtransect],size(SM(1).X1D,2),size(SM,2))';
% % 
% % TF=isoutlier(Z);
% % Zf=Z;Zf(TF == 1)=NaN;
% % Zmin=min(Zf);
% 
% figure('position',[136         156        1227         582]);
% plot(X2D(:),Z1D(:),'k.',X2Df(:),Zf(:),'r.')
% hold on;plot(X1D,Z1D(end,:),'g.');
% zmed=median(Zf,'omitnan');
% plot(X1Dt,zmed,'y-','linewidth',2);
% title(SA(end).File)
% plot([X0BeachOnly X0BeachOnly],[0 3],'k-');
% set(gca,'xdir','reverse')

% survey datetimes
Zdt=datetime([SM.Datenum],'convertfrom','datenum');

% year-month means
n=0;
for y=year(Zdt(1)):year(Zdt(end))
   for m=1:12
       n=n+1;
       Zyear(n)=y;
       Zmonth(n)=m;
       Zym(n,:)=Zf(1,:)*NaN;
       idx=find(year(Zdt) == y & month(Zdt) == m);
       if numel(idx) == 1
           Zym(n,:)=Zf(idx,:);
       elseif numel(idx) > 1
           Zym(n,:)=mean(Zf(idx,:),'omitnan');
       end
   end
end

% monthly means across all years
n=0;
for m=1:12
       n=n+1;
       Zm(n,:)=Zym(1,:)*NaN;
       idx=find(Zmonth == m);
       if numel(idx) == 1
           Zm(n,:)=Zym(idx,:);
       elseif numel(idx) > 1
           Zm(n,:)=mean(Zym(idx,:),'omitnan');
       end
end


% 3-month running means across all years
n=0;
for m3=1:12
       n=n+1;
       Zm3(n,:)=Zym(1,:)*NaN;
       if m3 == 1
           m=[1 2 12];
       elseif m3 == 12
           m=[1 11 12];
       else
           m=m3-1:m3+1;
       end

       idx=find(Zmonth == m(1) | Zmonth == m(2) | Zmonth == m(3));
       if numel(idx) == 1
           Zm3(n,:)=Zym(idx,:);
       elseif numel(idx) > 1
           Zm3(n,:)=mean(Zym(idx,:),'omitnan');
       end
end


% quarterly means from monthly means
n=0;
for q=1:4
       n=n+1;
       Zq(n,:)=Zm(1,:)*NaN;      
       idx=(q-1)*3+1:q*3;
       Zq(n,:)=mean(Zm(idx,:),'omitnan');
end

% figure;
% hold on;
% for m=1:4
%     plot(X1Dt,Zq(m,:),'-','linewidth',2)
%     hold on;
% end
% legend

% seasonal means from quarterly means
Zs(1,:)=mean(Zq(1:2,:),'omitnan');
Zs(2,:)=mean(Zq(3:4,:),'omitnan');

% figure;
% hold on;
% for m=1:2
%     plot(X1Dt,Zs(m,:),'-','linewidth',2)
%     hold on;
% end
% legend
    
% global mean
Zg=mean(Zs,'omitnan');


end