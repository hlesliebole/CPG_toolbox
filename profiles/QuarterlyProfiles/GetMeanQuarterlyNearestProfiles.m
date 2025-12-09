function [X0BeachOnly,X1Dt,QY,Q,Zyq,NS]=GetMeanQuarterlyNearestProfiles(SA,Syear)

%% returns year-quarter mean subaerial profiles for the input SA struct array
%
%  X0BeachOnly = Mop transect xshore location of truck back beach boundary
%  X1Dt(N) = N profile xshore x values (1m res) relative to Mop back beach point
%  QY(M) = M Quartery profile years
%  Q(M) = M Quartery profile quarters (1-4)
%  Zyq(M,N) = M QY(M)/Q(M) year/quarter mean profiles (Jan-Mar,Apr-Jun,Jul-Sep,Oct-Dec)
%  NS(M) = M QY(M)/Q(M) year/quarter number of surveys in quarter mean 

X0BeachOnly=NaN;
X1Dt=NaN;
QY=NaN;
Q=NaN;
Zyq=NaN;
NS=NaN;

%% get nearest point profile
Ytol=25; % 25m alongcoast nearest point tolerance
Xtol=5; % 5m cross shore gap interpolation tolerance
ZswashTol=0.05; % positive elevation change tolerance when stepping seaward
                %  on the profile below MSL, profile is cutoff if exceeded

%% convert survey dayes to their corresponding beach year
BeachYr=year(datetime([SA.Datenum]+92,'convertfrom','datenum'));

% reduce to only Truck LiDAR surveys in the desired beach survey year 
trk=find((strcmp({SA.Source},'Trk') | strcmp({SA.Source},'AtvMR') ) & ...
    BeachYr == Syear );
if ~isempty(trk)
SA=SA(trk);

% % Skip SfMdrone data for now
% sfm=find(strcmp({SA.Source},'SfMdrone'));
% if ~isempty(sfm)
%     SA(sfm)=[];
% end
% 
% % Skip AtvMR data for now
% amr=find(strcmp({SA.Source},'AtvMR'));
% if ~isempty(amr)
%     SA(amr)=[];
% end

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
Zf(Zf < -0.5)=NaN; % keep subaerial part
[zmin,imin]=min(Zf');
for n=1:size(Zf,1)
    if imin(n) > 0
        Zf(n,imin(n)+1:end)=NaN;
    end
end

[ibad,jbad]=find(Zf(:,1:end-2) < 0.774 & diff(Zf,2,2) > ZswashTol);
if ~isempty(ibad)
    for n=1:numel(ibad)
      Zf(ibad(n),jbad(n)+2:end)=NaN;
    end
end

% remove any additional outliers
TF=isoutlier(Zf,"mean");
Zf(TF == 1)=NaN;

% survey datetimes
Zdt=datetime([SA.Datenum],'convertfrom','datenum');

% year-month means
n=0;
for y=year(Zdt(1)):year(Zdt(end))
   for m=1:12
       n=n+1;
       Zyear(n)=y;
       Zmonth(n)=m;
       Zym(n,:)=Zf(1,:)*NaN;
       NSm(n,:)=Zf(1,:)*NaN;
       idx=find(year(Zdt) == y & month(Zdt) == m);
       if numel(idx) == 1
           Zym(n,:)=Zf(idx,:);
           NSm(n)=1;
       elseif numel(idx) > 1
           %Zym(n,:)=mean(Zf(idx,:),'omitnan');
           Zym(n,:)=mean(Zf(idx,:));
           %Zym(n,:)=median(Zf(idx,:),'omitnan');
           NSm(n)=numel(idx);
       end
   end
end

% reduce to year-quarter means by making a 3 month mean through the 
%   monthly profile and picking out the center month values of the
%   4 quarters
%Zym=movmedian(Zym,3,'omitnan');
Zym=movmean(Zym,3,'omitnan');
NSm=movsum(NSm,3,'omitnan');
idx=find(ismember(Zmonth,[2 5 8 11]));
Zyq=Zym(idx,:);
NS=NSm(idx);
Zyear=Zyear(idx);Zmonth=Zmonth(idx);
% reduce to year-quarter means with subaerial data
Zyq(Zyq < 0)=NaN;
idx=find(sum(Zyq,2,'omitnan') > 0);
Zyq=Zyq(idx,:);
QY=Zyear(idx);Q=Zmonth(idx);
% turn quarterly month number into a quarter number
Q=round((Q+1)/3);

end

end