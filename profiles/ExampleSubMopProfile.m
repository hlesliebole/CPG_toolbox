%SA=SAcombineMops(637,667);
SA=SAcombineMops(655,656);
SubMopNumber=655.0;
XgapTol=5;
YdistTol=10;

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans]=...
  GetCpgNearestPointSubMopProfiles(SA,SubMopNumber,XgapTol,YdistTol);

idx=find(contains({SA.File},'SolanaMixed'));

figure;plot(X1Dmop,Z1Dtrans(idx,:))