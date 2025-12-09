% creates mean quarterly subaerila beach profiles for SD
%  state beach mops

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
        535 607
        664 683
        684 706
        717 725
        740 757
        762 819
        825 854];

% loop thru state beaches

for sb=3%1:1

% loop through mop range for this beach
for MopNum=MopRng(sb,1)%:MopRng(sb,2)

% load SA mat file
load(['M' num2str(MopNum,'%5.5i') 'SA.mat'],'SA');

% get quarterly mean profiles
[X0BeachOnly,X1Dt,QY,Q,Zyq]=GetMeanQuarterlyNearestProfiles(SA);

for yq=1:numel(QY)
    ofile=[bname{sb} '_D' num2str(MopNum,'%4.4i') '_'...
        num2str(QY(yq)) '_' num2str(Q(yq)) '.dat'];
    fprintf('%s\n',ofile)
end

end
end