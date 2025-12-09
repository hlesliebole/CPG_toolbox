% make a combined SG grid for the malibu mop range
%  of L0891 to L0923
load MopTableUTM

Mop1=find(strcmp('L0791',table2cell(Mop(:,1))));
Mop2=find(strcmp('L0800',table2cell(Mop(:,1))));
MopM=find(strcmp('L0794',table2cell(Mop(:,1))));

% CG=SAcombineMops(Mop1,Mop2); % make combined grid struct array
% 
% [X,Y,Z]=SG2grid(CG,2);
figure('position',[ 211         169        1031         581]);%surf(X,Y,Z);shading flat;BeachColorbar;view(2)
hold on;
for n=Mop1:Mop2%CG(1).Mopnum'
    %PlotLabelMopTransectUTM(n,'3d','k','ShadeOff')
    PlotLabelLandMopTransect(n,'2d','m','ShadeOff')
end
PlotLabelLandMopTransect(MopM,'2d','y','ShadeOff')
plot_google_map('MapType', 'satellite');
set(gca,'xlim',[-118.5821 -118.5571]);
set(gca,'ylim',[ 34.0325   34.0448],'fontsize',16);
title('Ratner Beach CDIP MOnitoring and Prediction (MOP) Transects','fontsize',18)
makepng('RatnerMopOverview.png')

%title(CG(1).File)