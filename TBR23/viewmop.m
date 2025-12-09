function [p1]=viewmop(MopNumber)

matfile=['M' num2str(MopNumber,'%5.5i') 'SM.mat'];
load(matfile,'SM');

% % identify jumbos
% jumbo=find(contains({SM.File},'umbo'));
% % reduce local SM struct array to just jumbo surveys
% SM=SM(jumbo);

figure('position',[142         204        1096         506]);
hold on;xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88)');
title(['Mop ' num2str(MopNumber) ': ' num2str(size(SM,2)) ' Transect Profiles'],'fontsize',18);
for n=1:size(SM,2)
    p1=plot(SM(n).X1D,SM(n).Z1Dtransect,'-','DisplayName',...
        [datestr(SM(n).Datenum) ' ' SM(n).Source]);hold on
end
set(gca,'xdir','reverse','fontsize',12);box on;grid on;
legend('location','eastoutside')
end