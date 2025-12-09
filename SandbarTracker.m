
% load M00553GM.mat
% figure;plot(GM.Z1Dmean,'k-');hold on;plot(GM.Z1Dmedian,'g-')
% figure;plot(detrend(GM.Z1Dmean,'omitnan'));hold on;%plot(movmean(detrend(SM(208).Z1Dmean,'omitnan'),11))
MopNumber=675;%553;
load(['M' num2str(MopNumber,'%5.5i') 'SM.mat'],'SM')
ijumbo=find(contains({SM.File},'umbo')); % find jumbo survey
Zmean=mean(vertcat(SM(ijumbo).Z1Dtransect),'omitnan'); % global mean jumbo profile
X=SM(1).X1D; % xshore distance
Zsaq=Zmean(Zmean < 0 & Zmean > -8); % subaqueous section of global mean profile
Xsaq=SM(1).X1D(Zmean < 0 & Zmean > -8); % subaqueous xshore distance
Zsaqdt=detrend(Zsaq,'omitnan'); % detrended subaqueous global mean profile
Zsaqdtmm=movmean(Zsaqdt,11);
GZsaqTrend=Zsaq-Zsaqdt;

figure('position',[138   122   613   647]);
set(0, 'DefaultLineLineWidth', 2);
plot(X,Zmean,Xsaq,Zsaq,Xsaq,Zsaqdtmm,Xsaq,GZsaqTrend,'k:')
set(gca,'xdir','reverse','fontsize',12);grid on;
legend(['Global Mean Profile of (' num2str(numel(ijumbo)) ') Jumbos'],'Global Mean -8m to 0m Subaqueous Profile',...
    'Detrended Global Mean Subaqueous Profile','Global Mean Subaqueous Profile Trend',...
    'location','southeast','fontsize',12)
hold on;
[pks,ipks]=findpeaks(Zsaqdtmm);
for n=1:numel(pks)
plot([Xsaq(ipks) Xsaq(ipks)],[Zsaq(ipks) Zsaqdtmm(ipks)],...
    'm.:','markersize',20,'DisplayName',['Detrended Bar Peak: x=' num2str(round(Xsaq(ipks(n)))) 'm ; z=' num2str(Zsaq(ipks(n)),'%5.2fm')])
end

xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88)');
title(['Mop ' num2str(MopNumber) ' Global Mean Jumbo Profile'],'fontsize',16);

pngfile=['Mop' num2str(MopNumber) 'GlobalMeanJumboSandbar.png'];
makepng(pngfile)

% process last jumbo using global trend line
njumbo=find(ijumbo==24);
%njumbo=42; % Oct 2015
%njumbo=43; % Jan 2016
%njumbo=numel(ijumbo); % latest jumbo
Z=SM(ijumbo(njumbo)).Z1Dtransect;
ix=find(ismember(X,Xsaq));
Zsaq=Z(ix);
Zsaqdt=Zsaq-GZsaqTrend;
Zsaqdtmm=movmean(Zsaqdt,11);
[pks,ipks,w,p]=findpeaks(Zsaqdtmm);
if ~isempty(pks) % delete small peaks relative the max peak prominence p
    ipks(p < 0.25*max(p))=[];
    pks(p < 0.25*max(p))=[];
end

figure('position',[138   122   613   647]);
set(0, 'DefaultLineLineWidth', 2);
plot(X,Z,Xsaq,Zsaqdt,Xsaq,GZsaqTrend,'k:')
set(gca,'xdir','reverse','fontsize',12);grid on;
legend(['Jumbo Profile ' datestr(SM(ijumbo(njumbo)).Datenum)],...
    'Detrended Jumbo Profile','Global Mean Subaqueous Profile Trend',...
    'location','southeast','fontsize',12)
hold on;
for n=1:numel(pks)
plot([Xsaq(ipks(n)) Xsaq(ipks(n))],[Zsaq(ipks(n)) Zsaqdtmm(ipks(n))],...
    'm.:','markersize',20,'DisplayName',['Detrended Bar Peak: x=' num2str(round(Xsaq(ipks(n)))) 'm ; z=' num2str(Zsaq(ipks(n)),'%5.2fm')])
end

xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88)');
title(['Mop ' num2str(MopNumber) ' Jumbo ' datestr(SM(ijumbo(njumbo)).Datenum)],'fontsize',16);

pngfile=['Mop' num2str(MopNumber) 'Jumbo' datestr(SM(ijumbo(njumbo)).Datenum,'YYYYmmDD') 'Sandbar.png'];
makepng(pngfile)

% 
% figure;plot(detrend(Zmean(Zmean < 0 & Zmean > -8),'omitnan'));hold on;
% plot(movmean(detrend(Zmean(Zmean < 0 & Zmean > -8)),11));
% plot(t,x,t,y,t,x-y,':k')
% legend('Input Data','Detrended Data','Trend','Location','northwest') 
% 
% figure;plot(detrend(SM(208).Z1Dmean,'omitnan'));hold on;plot(movmean(detrend(SM(208).Z1Dmean,'omitnan'),11))