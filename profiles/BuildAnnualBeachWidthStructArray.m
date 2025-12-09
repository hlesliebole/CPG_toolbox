% creates mean quarterly subaerila beach profiles for SD
%  state beach mops
close all
clearvars

% Syear=2023; % Survey year to reduce to quarterly means
% % open mop prolile definition file
% fid=fopen('MopLatLonEutnNutmOrientation.dat','w');

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
bname{10}='ImperialBeach';
bname{11}='DelMar';
bname{12}='SolanaBeach';
bname{13}='Oceanside';

% define 9 mop ranges for the state beaches
MopRng=[2 28 
        85 157
        536 607
        664 683
        684 706
        717 725
        740 757
        762 819
        825 854
        39 61
        608 635
        638 663
        866 903];

% loop thru  beaches

m=0;
for sb=1:13

% loop through mop range for this beach
for MopNum=MopRng(sb,1):MopRng(sb,2)
%%
fprintf('%i %i\n',sb,MopNum)

% load SA mat file
load(['M' num2str(MopNum,'%5.5i') 'SA.mat'],'SA');

NumSubTrans=5; % 1 will just return the main transect info; 100 would make
                % statistical profiles with complete use of any 1m res LiDAR
XgapTol=5; % patch any cross-shore profile gaps < 5m wide
YdistTol=25; % can go as large as 50m up- downcoast for a single Mop area

% get nearest point profiles

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);

% get annual beach widths
[TZdatetime,TZ1Dmon,TZ1Dqtr,TZ1Dsea,TZ1Dann,TZ1Dglo]=...
                     TimeAveragedMeanProfiles(Zdatetime,Z1Dmedian);
%%
Znavd88=1.344;%Znavd88=0.774;
[BWmin,BWmax]=GetCpgProfileBeachWidths(Znavd88,X1Dmop,TZ1Dann);

m=m+1;
ABW(m).BeachNum=sb;
ABW(m).Mopnum=MopNum;
ABW(m).Byears=TZdatetime.ann;
ABW(m).BWmin=BWmin;
ABW(m).BWmax=BWmax;

%figure;plot(TZdatetime.ann,BWmax,'.-')
%%

%% get quarterly mean profiles

%  X0BeachOnly = Mop xshore location of truck back beach boundary
%  X1Dt(N) = N profile xshore x values (1m res) relative to Mop back beach point
%  QY(M) = M Quartery profile years
%  Q(M) = M Quartery profile quarters (1-4)
%  Zyq(M,N) = M QY(M)/Q(M) year/quarter mean profiles (Jan-Mar,Apr-Jun,Jul-Sep,Oct-Dec)
%  NS(M) = M QY(M)/Q(M) year/quarter number of surveys in quarter mean 
% 
% [X0BeachOnly,X1Dt,QY,Q,Zyq,NS]=GetMeanQuarterlyNearestProfiles(SA,Syear);
% 
% %  Get utm and latlon coords of xshore profile points
% [Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(MopNum,X1Dt);
% 
% %figure;
% for yq=1:numel(QY)
%     %plot3(X1Dt,(QY(yq)+Q(yq)/4)*ones(size(X1Dt)),Zyq(yq,:),'-');hold on;
%     ofile=[bname{sb} '/' bname{sb} '_D' num2str(MopNum,'%4.4i') '_'...
%         num2str(QY(yq)) '_' num2str(Q(yq)) '.dat'];
%     fprintf('%s\n',ofile)
%     % find last valid data point in the profile
%     ilast=find(~isnan(Zyq(yq,:)),1,'last');
%     fid=fopen(ofile,'w');
%     for m=1:ilast
%         %fprintf(fid,'%i %8.3f\n',n-1,Zyq(yq,n));
%         %for m=1:numel(X1Dt)
%         fprintf(fid,'%i %10.1f %10.1f %13.7f %13.7f %8.2f\n',...
%             X1Dt(m),Xutm(m),Yutm(m),Lon(m),Lat(m),Zyq(yq,m));
%         %end
%     end
%     fclose(fid);
% end

end
end

%save AnnualBeachWidths.mat ABW