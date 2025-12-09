%function ProfileXshoreGapAssessment(MopNumber)
clearvars
MopNumber=674;

%load the SA struct array
load(['M' num2str(MopNumber,'%5.5i') 'SA.mat'],'SA')

% first find all the jumbo of jetski only survey data indexes
jdx=find(contains({SA.File},'umbo') | contains({SA.File},'etski'));

% reduce SA struct array to just these surveys
SA=SA(jdx);

% now pass to the profile function using a modest YdistTol that should
%  pull in most jumbo data that is near the transect
NumSubTrans=1; % 1 will just return the main transect info; 100 would make
                % statistical profiles with complete use of any 1m res LiDAR
XgapTol=0; % don't fill gaps for gap assessment purposes
YdistTol=20; % allow some tolerance for being off the Mop line

% get nearest point profiles
[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);
% Z1Dtrans is the matrix

% also want the temporal long-term global median profile for gap-depth assessment
[TZdatetime,TZ1Dmon,TZ1Dqtr,TZ1Dsea,TZ1Dann,TZ1Dglo]=...
                     TimeAveragedMeanProfiles(Zdatetime,Z1Dmedian);

%  loop through all the Z1Dtrans profiles and collect their gapsize
%    information as a function of the gap's global median depth

%  Use third party matlab code for this

%figure;plot(X1Dcpg,TZ1Dglo,'b-')

for n=1:size(Z1Dtrans,1)
    gs=gapsize(Z1Dtrans(n,:)); % gap value for every cross shore point
    gs=diff(gs); % isolate single gap value for each gap with diff
    idx=find(gs > 0); % find shoreward side x of each gap
    ndx=idx+round(gs(idx)/2); % estimate x of middle of gap   
    if ndx(end) > size(Z1Dtrans,2) % avoid offshore overshoot of last gap
        ndx(end) = size(Z1Dtrans,2);
    end 
    zgs=TZ1Dglo(ndx); % assign a depth associated with the gap
    % save in struct array
    GAP(n).Datenum=SA(n).Datenum;
    GAP(n).GapSize=gs(idx);
    GAP(n).GapDepth=zgs;
end

% plot histograms
figure;histogram([GAP.GapSize],.5:1:30.5);
figure;plot([GAP.GapDepth],[GAP.GapSize],'*');set(gca,'ylim',[0 20])

figure;histogram2([GAP.GapDepth],[GAP.GapSize],-10:1:5,.5:1:30.5,'DisplayStyle','tile');
colormap(jet(64))
%%
figure('position',[78          72        1138         702]);
histogram2([GAP.GapDepth],[GAP.GapSize],-10.5:1:4.5,2.5:1:30.5,'FaceColor','flat','facealpha',.8);
colormap(jet(64))
set(gca,'view',[170.3431   35.2294],'xtick',-10:1:5,'ylim',[1 30],'ytick',2:2:30)
set(gca,'fontsize',18);xlabel('Mean Elevation (m,NAVD88)');ylabel('Gap Size (m)');zlabel('N');
colormap(jet(64));colorbar
title('Gap Size >= 2m Vs. Mean Elevation')
set(gca,'linewidth',2)


function sz=gapsize(x)
% Calculates the gap size of a vector. 
%
% A gap is defined as the number of consequtive nans.
%
% USAGE: sz=gapsize(x)
% sz is same length as input. 
%
% example:
% x=rand(20,1);
% x(x>.5)=nan; 
% [x gapsize(x)]
%
% aslak grinsted 2010
x=~isnan(x);
hasdata=[0;find(x(:)); length(x)+1];
sz=zeros(size(x));
for ii=1:length(hasdata)-1
    ix=hasdata(ii)+1:hasdata(ii+1)-1;
    sz(ix)=length(ix);
end

end







