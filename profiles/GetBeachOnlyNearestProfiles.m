function [X0BeachOnly,X1Dbo,Zbo,Zavgbo,Zqbo,Zsbo,Zgbo]=GetBeachOnlyNearestProfiles(SA,Ytol,Xtol)

%% returns mean profiles for the input SA struct array
%
%  X0BeachOnly = Mop xshore location of truck back beach boundary
%  X1Dbo(N) = N xshore x values (1m res) relative to X0BeachOnly
%  Zbo(M,N) = M survey profiles relative to X0BeachOnly
%  Zavgbo(12,N) = monthly mean profiles (calculated by year first and then across all years)
%  Zqbo(4,N) = quarterly mean profiles (Jan-Mar,Apr-Jun,Jul-Sep,Oct-Dec)
%  Zsbo(2,N) = seasonal mean profiles (Jan-Jun,Jul-Dec)
%  Zgbo(1,N) = global mean profile

%% get nearest point profile
% Ytol=50; % 50m alongcoast nearest point tolerance
% Xtol=10; % 5m cross shore gap interpolation tolerance

[X1Dbo,Z1D]=GetAllNearestPointsProfiles(SA,Ytol,Xtol);
% size(X1Dbo)
% size(Z1D)

%X2D=repmat(X1Dbo,size(Z1D,1),1);
%% find the truck back beach boundary along the mop tranect
trk=find(strcmp({SA.Source},'Trk'));
if ~isempty(trk)
for n=1:numel(trk)
    xback(n)=X1Dbo(find(~isnan(Z1D(trk(n),:)), 1 ));
end
else
    xback=0;
end
X0BeachOnly=median(xback);
% fprintf('Truck Back Beach x= %6.1f \n',X0BeachOnly)

%% reduce profiles to common back beach x location based on truck beach-only data
idx=find(X1Dbo >= X0BeachOnly);
Z1Dt=Z1D(:,idx);
X1Dbot=X1Dbo(idx);

%% remove any elevations seaward of the min profile elevation
%  (lidar swash filter)
 Zbo=Z1Dt;
[ Zavgboin,imin]=min( Zbo');
for n=1:size( Zbo,1)
    if imin(n) > 0
         Zbo(n,imin(n)+1:end)=NaN;
    end
end

% remove any additional outliers
TF=isoutlier( Zbo,"mean");
 Zbo(TF == 1)=NaN;

%X2Df=repmat(X1Dbot,size( Zbo,1),1);

% 
% 
% Z= Zbo;
% Z(Z > zc)=NaN;
% [ Zavgboax,imax]=max(Z');
% imax(isnan( Zavgboax))=NaN;
% imax(zc- Zavgboax > 0.1)=NaN;
%  Zavgboax(zc- Zavgboax > 0.1)=NaN;
% figure;
% subplot(2,1,1);plot( Zavgboax);hold on;subplot(2,1,2);plot(imax);hold on;plot(imax,'o')
% 
% Z= Zbo;
% Z(Z < zc)=NaN;
% [ Zavgboin,imin]=min(Z');
% imin(isnan( Zavgboin))=NaN;
% imin( Zavgboin-zc > 0.1)=NaN;
%  Zavgboin( Zavgboin-zc > 0.1)=NaN;
% %figure;
% subplot(2,1,1);plot( Zavgboin);subplot(2,1,2);plot(imin);hold on;plot(imin,'^')
% imean=(nanmean(vertcat(imin,imax)));plot(imean,'+')
% fprintf('%i surveys have valid %7.3f m contour location data.\n',numel(find(~isnan(imean))),zc)
% dztol=0.1;
% Z= Zbo;
% %Z(Z < zc)=NaN;
% [ Zavgboin,imin]=min(abs(Z-zc)');
% imin(isnan( Zavgboin))=NaN;
% imin( Zavgboin > dztol)=NaN;
% % Zavgboin( Zavgboin-zc > 0.1)=NaN;
% figure
% subplot(2,1,1);plot( Zavgboin);hold on;plot([0 numel( Zavgboin)],[dztol dztol],'k--');
% set(gca,'yscale','log');
% xlabel('Survey Number');ylabel('Dist from contour (m)');
% title(['Min Profile Point Distance from ' num2str(zc) ' m contour']);
% subplot(2,1,2);plot(imin(~isnan(imin)));hold on;plot(X1Dbot(imin(~isnan(imin))),'s')
% fprintf('%i surveys have a valid %6.3f m contour location data.\n',numel(find(~isnan(imin))),zc)
% xlabel('Survey Number');ylabel('Dist from back beach (m)');
% title(['Contour Distance from (Truck) back beach (m)' num2str(zc) ' m contour']);
% % load M00582SM.mat
% % Z=reshape([SM.Z1Dtransect],size(SM(1).X1Dbo,2),size(SM,2))';
% % 
% % TF=isoutlier(Z);
% %  Zbo=Z; Zbo(TF == 1)=NaN;
% %  Zavgboin=min( Zbo);
% 
% figure('position',[136         156        1227         582]);
% plot(X2D(:),Z1D(:),'k.',X2Df(:), Zbo(:),'r.')
% hold on;plot(X1Dbo,Z1D(end,:),'g.');
%  Zavgboed=median( Zbo,'omitnan');
% plot(X1Dbot, Zavgboed,'y-','linewidth',2);
% title(SA(end).File)
% plot([X0BeachOnly X0BeachOnly],[0 3],'k-');
% set(gca,'xdir','reverse')

% survey datetimes
Zdt=datetime([SA.Datenum],'convertfrom','datenum');

% year-month means
n=0;
for y=year(Zdt(1)):year(Zdt(end))
   for m=1:12
       n=n+1;
       Zyear(n)=y;
        Zavgboonth(n)=m;
       Zym(n,:)= Zbo(1,:)*NaN;
       idx=find(year(Zdt) == y & month(Zdt) == m);
       if numel(idx) == 1
           Zym(n,:)= Zbo(idx,:);
       elseif numel(idx) > 1
           Zym(n,:)=mean( Zbo(idx,:),'omitnan');
       end
   end
end

% monthly means across all years
n=0;
for m=1:12
       n=n+1;
        Zavgbo(n,:)=Zym(1,:)*NaN;
       idx=find( Zavgboonth == m);
       if numel(idx) == 1
            Zavgbo(n,:)=Zym(idx,:);
       elseif numel(idx) > 1
            Zavgbo(n,:)=mean(Zym(idx,:),'omitnan');
       end
end


% 3-month running means across all years
n=0;
for m3=1:12
       n=n+1;
        Zavgbo3(n,:)=Zym(1,:)*NaN;
       if m3 == 1
           m=[1 2 12];
       elseif m3 == 12
           m=[1 11 12];
       else
           m=m3-1:m3+1;
       end

       idx=find( Zavgboonth == m(1) |  Zavgboonth == m(2) |  Zavgboonth == m(3));
       if numel(idx) == 1
            Zavgbo3(n,:)=Zym(idx,:);
       elseif numel(idx) > 1
            Zavgbo3(n,:)=mean(Zym(idx,:),'omitnan');
       end
end


% quarterly means from monthly means
n=0;
for q=1:4
       n=n+1;
       Zqbo(n,:)= Zavgbo(1,:)*NaN;      
       idx=(q-1)*3+1:q*3;
       Zqbo(n,:)=mean( Zavgbo(idx,:),'omitnan');
end

% figure;
% hold on;
% for m=1:4
%     plot(X1Dbot,Zqbo(m,:),'-','linewidth',2)
%     hold on;
% end
% legend

% seasonal means from quarterly means
Zsbo(1,:)=mean(Zqbo(1:2,:),'omitnan');
Zsbo(2,:)=mean(Zqbo(3:4,:),'omitnan');

% figure;
% hold on;
% for m=1:2
%     plot(X1Dbot,Zsbo(m,:),'-','linewidth',2)
%     hold on;
% end
% legend
    
% global mean
Zgbo=mean(Zsbo,'omitnan');


end