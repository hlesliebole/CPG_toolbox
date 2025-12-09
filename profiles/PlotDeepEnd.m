
figure;
load M00582SM.mat
[X0BeachOnly,X1Dt,Zf,Zm,Zm3,Zq,Zs,Zg]=GetMeanMedianProfiles(SM);

plot([SM.Datenum],Zf(:,400),'+');hold on;plot([SM.Datenum],Zf(:,450),'+');hold on;plot([SM.Datenum],Zf(:,500),'+');plot([SM.Datenum],Zf(:,550),'+')
datetick
grid on;

load M00616SM.mat
%load M00672SM.mat
[X0BeachOnly,X1Dt,Zf,Zm,Zm3,Zq,Zs,Zg]=GetMeanMedianProfiles(SM);

% plot([SM.Datenum],Zf(:,400),'o');hold on;plot([SM.Datenum],Zf(:,450),'o');hold on;plot([SM.Datenum],Zf(:,500),'o');plot([SM.Datenum],Zf(:,550),'o')
% datetick

plot([SM.Datenum],Zf(:,500),'o');hold on;plot([SM.Datenum],Zf(:,550),'o');hold on;plot([SM.Datenum],Zf(:,600),'o');plot([SM.Datenum],Zf(:,650),'o')
datetick

