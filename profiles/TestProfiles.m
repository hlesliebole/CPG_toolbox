% test of profile code
close all
clearvars

NumSubTrans=11;
XgapTol=15;
YdistTol=25;

load M00582SA.mat

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);

figure;
%imagesc(X1Dcpg,datenum(Zdatetime),Z1Dmean,'alphadata',~isnan(Z1Dtrans));
%datetick('y');
pcolor(X1Dcpg,Zdatetime,Z1Dmean);shading flat;
BeachColorbar;