% make a combined SG grid for the malibu mop range
%  of L0891 to L0923

Mop1=find(strcmp('L0891',table2cell(Mop(:,1))));
Mop2=find(strcmp('L0923',table2cell(Mop(:,1))));

CG=SGcombineMops(Mop1,Mop2); % make combined grid struct array

[X,Y,Z]=SG2grid(CG,1);
figure;surf(X,Y,Z);shading flat;BeachColorbar;view(2)
hold on;
for n=CG(1).Mopnum'
    PlotLabelMopTransectUTM(n,'3d','k','ShadeOff')
end
title(CG(1).File)