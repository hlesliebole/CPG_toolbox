% creates mean quarterly subaerila beach profiles for SD
%  state beach mops
close all
clearvars

Syear=2024; % Survey year to reduce to quarterly means
% open mop prolile definition file
fid=fopen('MopLatLonEutnNutmOrientation.dat','w');

% define 9 state beach names for files
bname{1}='BorderField';
bname{2}='SilverStrand';
bname{3}='TorreyPines';
bname{4}='Cardiff';
bname{5}='SanElijo';
bname{6}='Moonlight';
bname{7}='Leucadia';
bname{8}='SouthCarlsbad';
bname{9}='Carlsbad';

% define 9 mop ranges for the state beaches
MopRng=[2 28 
        85 157
        536 607
        664 683
        684 706
        717 725
        740 757
        762 819
        825 854];

% loop thru state beaches

for sb=1:9

% loop through mop range for this beach
for MopNum=MopRng(sb,1):MopRng(sb,2)

% load SA mat file
load(['M' num2str(MopNum,'%5.5i') 'SA.mat'],'SA');

%% get quarterly mean profiles

%  X0BeachOnly = Mop xshore location of truck back beach boundary
%  X1Dt(N) = N profile xshore x values (1m res) relative to Mop back beach point
%  QY(M) = M Quartery profile years
%  Q(M) = M Quartery profile quarters (1-4)
%  Zyq(M,N) = M QY(M)/Q(M) year/quarter mean profiles (Jan-Mar,Apr-Jun,Jul-Sep,Oct-Dec)
%  NS(M) = M QY(M)/Q(M) year/quarter number of surveys in quarter mean 

[X0BeachOnly,X1Dt,QY,Q,Zyq,NS]=GetMeanQuarterlyNearestProfiles(SA,Syear);

%  Get utm and latlon coords of xshore profile points
[Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(MopNum,X1Dt);

%figure;
for yq=1:numel(QY)
    %plot3(X1Dt,(QY(yq)+Q(yq)/4)*ones(size(X1Dt)),Zyq(yq,:),'-');hold on;
    ofile=[bname{sb} '/' bname{sb} '_D' num2str(MopNum,'%4.4i') '_'...
        num2str(QY(yq)) '_' num2str(Q(yq)) '.dat'];
    fprintf('%s\n',ofile)
    % find last valid data point in the profile
    ilast=find(~isnan(Zyq(yq,:)),1,'last');
    fid=fopen(ofile,'w');
    for m=1:ilast
        %fprintf(fid,'%i %8.3f\n',n-1,Zyq(yq,n));
        %for m=1:numel(X1Dt)
        fprintf(fid,'%i %10.1f %10.1f %13.7f %13.7f %8.2f\n',...
            X1Dt(m),Xutm(m),Yutm(m),Lon(m),Lat(m),Zyq(yq,m));
        %end
    end
    fclose(fid);
end

end
end