%function [X0BeachOnly,X1Dt,QY,Q,Zyq]=GetMeanQuarterlyNearestProfiles(SA)

%% returns year-quarter mean subaerial profiles for the input SA struct array
%
%  X0BeachOnly = Mop xshore location of truck back beach boundary
%  X1Dt(N) = N profile xshore x values (1m res) relative to X0BeachOnly
%  QY(M) = M Quartery profile years
%  Q(M) = M Quartery profile quarters (1-4)
%  Zyq(M,N) = M QY(M)/Q(M) year/quarter mean profiles (Jan-Mar,Apr-Jun,Jul-Sep,Oct-Dec)


%% get nearest point profile
Ytol=25; % 25m alongcoast nearest point tolerance
Xtol=5; % 5m cross shore gap interpolation tolerance

% Skip SfMdrone data for now
sfm=find(strcmp({SA.Source},'SfMdrone'));
if ~isempty(sfm)
    SA(sfm)=[];
end

% Skip AtvMR data for now
amr=find(strcmp({SA.Source},'AtvMR'));
if ~isempty(amr)
    SA(amr)=[];
end

[X1D,Z1D]=GetNearestPointsProfiles(SA,Ytol,Xtol);

%X2D=repmat(X1D,size(Z1D,1),1);
%% find the truck back beach boundary along the mop tranect
trk=find(strcmp({SA.Source},'Trk'));
if ~isempty(trk)
for n=1:numel(trk)
    idx=find(~isnan(Z1D(trk(n),:)), 1 );
    if ~isempty(idx)
        xback(n)=X1D(find(~isnan(Z1D(trk(n),:)), 1 ));
    else
        xback(n)=NaN;
    end
end
else
    xback=0;
end
X0BeachOnly=median(xback,'omitnan');
if isnan(X0BeachOnly)
    X0BeachOnly=0;
end

% fprintf('Truck Back Beach x= %6.1f \n',X0BeachOnly)

%% reduce profiles to common back beach x location based on truck beach-only data
idx=find(X1D >= X0BeachOnly);
Z1Dt=Z1D(:,idx);
X1Dt=X1D(idx);

%% remove any elevations seaward of the global min profile elevation
%  (lidar swash filter)
Zf=Z1Dt;
Zf(Zf < 0)=NaN; % keep subaerial part
[zmin,imin]=min(Zf');
for n=1:size(Zf,1)
    if imin(n) > 0
        Zf(n,imin(n)+1:end)=NaN;
    end
end

% find vertical max and seaward x max values in data set
zmax=max(Zf(:));
[ix,iy]=find(~isnan(Zf));
xmax=max(iy);

% remove any additional outliers
TF=isoutlier(Zf,"mean");
Zf(TF == 1)=NaN;

% survey datetimes
Zdt=datetime([SA.Datenum],'convertfrom','datenum');

% year-month means
figure
n=0;
for y=year(Zdt(1)):year(Zdt(end))
   for m=1:12
       n=n+1;
       Zyear(n)=y;
       Zmonth(n)=m;
       Zym(n,:)=Zf(1,:)*NaN;
       idx=find(year(Zdt) == y & month(Zdt) == m);
       if ~isempty(idx)
           for i=idx
               plot(X1Dt,Zf(i,:),'-');hold on;
           end
           title([num2str(y) num2str(m)])
           set(gca,'xlim',[X1Dt(1) X1Dt(xmax)],'ylim',[0 zmax]);grid on;
           pause
           clf
       end

       if numel(idx) == 1
           Zym(n,:)=Zf(idx,:);
       elseif numel(idx) > 1
           %Zym(n,:)=mean(Zf(idx,:),'omitnan');
           Zym(n,:)=median(Zf(idx,:),'omitnan');
       end
   end
end

% reduce to year-quarter means
Zym=movmedian(Zym,3,'omitnan');
%Zym=movmean(Zym,3,'omitnan');
idx=find(ismember(Zmonth,[2 5 8 11]));
Zyq=Zym(idx,:);
Zyear=Zyear(idx);Zmonth=Zmonth(idx);
% reduce to year-quarter means with subaerial data
Zyq(Zyq < 0)=NaN;
idx=find(sum(Zyq,2,'omitnan') > 0);
Zyq=Zyq(idx,:);
QY=Zyear(idx);Q=Zmonth(idx);
% turn quesrterly month number into a quarter number
Q=round((Q+1)/3);

%end