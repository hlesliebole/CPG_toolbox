
MopNumber=550;%553;%675;%%553;582;
load(['M' num2str(MopNumber,'%5.5i') 'SM.mat'],'SM')
idrone=find(contains({SM.Source},'RTKdrone')); % find RTK drone
SM(idrone)=[]; % remove RRTdrone for now
idrone=find(contains({SM.Source},'AtvMR')); % find atv miniranger
SM(idrone)=[]; % remove AtvMR for now
ijumbo=find(contains({SM.File},'umbo')); % find jumbo survey
Zmean=mean(vertcat(SM(ijumbo).Z1Dtransect),'omitnan'); % global mean jumbo profile
Zmean=mean(vertcat(SM.Z1Dtransect),'omitnan'); % global mean jumbo profile
X=SM(1).X1D; % xshore distance
Zmean(X < -20)=NaN;

% subaerial profile
Zsa=Zmean(Zmean >= 0 ); % subaerial section of global mean profile
Xsa=SM(1).X1D(Zmean >= 0); % subaerial xshore distance
Zsadt=detrend(Zsa,'omitnan'); % detrended subaerial global mean profile
Zsadtmm=movmean(Zsadt,11);
GZsaTrend=Zsa-Zsadt;

% subaqueous profile
Zsaq=Zmean(Zmean < 0 & Zmean > -8); % subaqueous section of global mean profile
Xsaq=SM(1).X1D(Zmean < 0 & Zmean > -8); % subaqueous xshore distance
Zsaqdt=detrend(Zsaq,'omitnan'); % detrended subaqueous global mean profile
Zsaqdtmm=movmean(Zsaqdt,11);
GZsaqTrend=Zsaq-Zsaqdt;

figure('position',[138   122   613   647]);
set(0, 'DefaultLineLineWidth', 2);
plot(X,Zmean,Xsaq,Zsaq,Xsaq,Zsaqdtmm,Xsa,Zsadtmm,Xsaq,GZsaqTrend,'k:',Xsa,GZsaTrend,'k--')
set(gca,'xdir','reverse','fontsize',12);grid on;
legend(['Global Mean Profile of (' num2str(numel(ijumbo)) ') Jumbos'],'Global Mean -8m to 0m Subaqueous Profile',...
    'Detrended Global Mean Subaqueous Profile','Detrended Global Mean Subaerial Profile',...
    'Global Mean Subaqueous Profile Trend',...
    'Global Mean Subaerial Profile Trend',...
    'location','southeast','fontsize',12)
hold on;
[pks,ipks,w,p]=findpeaks(Zsaqdtmm);
if ~isempty(pks) % delete small peaks relative the max peak prominence p
    ipks(p < 0.25*max(p))=[];
    pks(p < 0.25*max(p))=[];
end
for n=1:numel(pks)
plot([Xsaq(ipks) Xsaq(ipks)],[Zsaq(ipks) Zsaqdtmm(ipks)],...
    'm.:','markersize',20,'DisplayName',['Detrended Bar Peak: x=' num2str(round(Xsaq(ipks(n)))) 'm ; z=' num2str(Zsaq(ipks(n)),'%5.2fm')])
end

[bpks,ibpks,w,p]=findpeaks(Zsadtmm);
if ~isempty(bpks) % delete small peaks relative the max peak prominence p
    ibpks(p < 0.25*max(p))=[];
    bpks(p < 0.25*max(p))=[];
end
for n=1:numel(bpks)
plot([Xsa(ibpks(n)) Xsa(ibpks(n))],[Zsa(ibpks(n)) Zsadtmm(ibpks(n))],...
    'g.:','markersize',20,'DisplayName',['Detrended Berm Peak: x=' num2str(round(Xsa(ibpks(n)))) 'm ; z=' num2str(Zsa(ibpks(n)),'%5.2fm')])
end


xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88)');
title(['Mop ' num2str(MopNumber) ' Global Mean Jumbo Profile'],'fontsize',16);

% pngfile=['Mop' num2str(MopNumber) 'GlobalMeanJumboSandbar.png'];
% makepng(pngfile)
% 

figure('position',[154 272 1075 634]);hold on;
for njumbo=1:numel(ijumbo)
    
Z=SM(ijumbo(njumbo)).Z1Dtransect;
ix=find(ismember(X,Xsaq));
Zsaq=Z(ix);
Zsaqdt=Zsaq-GZsaqTrend;
Zsaqdtmm=movmean(Zsaqdt,11);
[pks,ipks,w,p]=findpeaks(Zsaqdtmm);
%[sp,npks]=sort(p,'descend');ipks=ipks(npks);pks=pks(npks);p=p(npks);
if ~isempty(pks) % delete small peaks relative the max peak prominence p
    ipks(p < 0.10*max(p))=[];
    pks(p < 0.10*max(p))=[];
else
    if min(SM(ijumbo(njumbo)).Z1Dtransect) < -4
    fprintf('No Bar %i %s\n',njumbo,datestr(SM(ijumbo(njumbo)).Datenum))
    end
end

sdate=day(datetime(SM(ijumbo(njumbo)).Datenum,'convertfrom','datenum'),'dayofyear');
col=['k','r','g'];

%subplot(3,1,1);hold on;plot([1 365],[min(Xsaq) min(Xsaq)],'k--')
Gmsl=max(Xsa(Xsa > 0.774));
for n=1:numel(pks)
    nn=n;if n > 3;nn=3;end
    subplot(3,1,1);hold on;datetick;grid on;
    plot(sdate,Xsaq(ipks(n))-Gmsl,'.','color',col(nn),'markersize',10);
    subplot(3,1,2);hold on;datetick;grid on;
    plot(sdate,Zsaq(ipks(n)),'.','color',col(nn),'markersize',10);
    subplot(3,1,3);hold on;datetick;grid on;
    plot(sdate,Zsaqdtmm(ipks(n)),'.','color',col(nn),'markersize',10)
end

end

%  Berm Plot
figure('position',[154 272 1075 634]);hold on;
%for njumbo=1:numel(ijumbo)
for ns=1:size(SM,2)
    
Z=SM(ns).Z1Dtransect;
ix=find(ismember(X,Xsa));
Zsa=Z(ix);
Zsadt=Zsa-GZsaTrend;
Zsadtmm=movmean(Zsadt,3);
[pks,ipks,w,p]=findpeaks(Zsadtmm);
%[sp,npks]=sort(p,'descend');ipks=ipks(npks);pks=pks(npks);p=p(npks);
if ~isempty(pks) % delete small peaks relative the max peak prominence p
%     ipks(p < 0.25*max(p))=[];
%     pks(p < 0.25*max(p))=[];
else
    [zmax,imax]=max(Zsa);
    fprintf('No Berm %i %s %5.1f %5.1f\n',ns,datestr(SM(ns).Datenum),X(ix(imax)),zmax)
    
end

sdate=day(datetime(SM(ns).Datenum,'convertfrom','datenum'),'dayofyear');
col=['k','r','g'];

%subplot(3,1,1);hold on;plot([1 365],[min(Xsaq) min(Xsaq)],'k--')
Gmsl=max(Xsa(Xsa > 0.774));
for n=1:numel(pks)
    nn=n;if n > 3;nn=3;end
    subplot(3,1,1);hold on;datetick;grid on;
    plot(sdate,Xsa(ipks(n))-Gmsl,'.','color',col(nn),'markersize',10);
    subplot(3,1,2);hold on;datetick;grid on;
    plot(sdate,Zsa(ipks(n)),'.','color',col(nn),'markersize',10);
    subplot(3,1,3);hold on;datetick;grid on;
    plot(sdate,Zsadtmm(ipks(n)),'.','color',col(nn),'markersize',10)
end

end



% % process last jumbo using global trend line
% njumbo=42; % Oct 2015
% %njumbo=43; % Jan 2016
% %njumbo=numel(ijumbo); % latest jumbo
% Z=SM(ijumbo(njumbo)).Z1Dtransect;
% ix=find(ismember(X,Xsaq));
% Zsaq=Z(ix);
% Zsaqdt=Zsaq-GZsaqTrend;
% Zsaqdtmm=movmean(Zsaqdt,11);
% [pks,ipks,w,p]=findpeaks(Zsaqdtmm);
% if ~isempty(pks) % delete small peaks relative the max peak prominence p
%     ipks(p < 0.25*max(p))=[];
%     pks(p < 0.25*max(p))=[];
% end
% 
% figure('position',[138   122   613   647]);
% set(0, 'DefaultLineLineWidth', 2);
% plot(X,Z,Xsaq,Zsaqdt,Xsaq,GZsaqTrend,'k:')
% set(gca,'xdir','reverse','fontsize',12);grid on;
% legend(['Jumbo Profile ' datestr(SM(ijumbo(njumbo)).Datenum)],...
%     'Detrended Jumbo Profile','Global Mean Subaqueous Profile Trend',...
%     'location','southeast','fontsize',12)
% hold on;
% for n=1:numel(pks)
% plot([Xsaq(ipks(n)) Xsaq(ipks(n))],[Zsaq(ipks(n)) Zsaqdtmm(ipks(n))],...
%     'm.:','markersize',20,'DisplayName',['Detrended Bar Peak: x=' num2str(round(Xsaq(ipks(n)))) 'm ; z=' num2str(Zsaq(ipks(n)),'%5.2fm')])
% end
% 
% xlabel('Xshore Distance (m)');ylabel('Elevation (m, NAVD88)');
% title(['Mop ' num2str(MopNumber) ' Jumbo ' datestr(SM(ijumbo(njumbo)).Datenum)],'fontsize',16);
% 
