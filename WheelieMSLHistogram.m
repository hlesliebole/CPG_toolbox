
load M00581SM.mat
SS=SMtoShoreline(581,581)
% 
% idx=find(strcmpi({SS.Source},'iG8wheel') == 1);
% rmsl=round([SS(idx).MSLbeachWidths]);
% figure
% histogram(rmsl,[25.5:60.5]);
% 
% idx=find(strcmpi({SS.Source},'iG8wheel') == 1 &...
% round([SS.MSLbeachWidths]) == 31);
% figure
% for n=idx
%     plot(-SM(n).X1D,SM(n).Z1Dmean)
%     hold on;
% end
    
idx=find(strcmpi({SS.Source},'iG8wheel') == 1);
rmsl=round([SS(idx).MHWbeachWidths]);
figure
histogram(rmsl);

idx=find(strcmpi({SS.Source},'iG8wheel') == 1 &...
round([SS.MHWbeachWidths]) == 16);
figure
for n=idx
    plot(-SM(n).X1D,SM(n).Z1Dmean)
    hold on;
end

idx=find(strcmpi({SS.Source},'iG8wheel') == 1);
v=[SS(idx).Volume];
figure
histogram(round(v));

idx=find(strcmpi({SS.Source},'iG8wheel') == 1 &...
round([SS.Volume]) == 19);
figure
xs=[];
ss=[];
for n=idx
    [Volume,Centroid,Slope,SlopeLF,rmseLF]=...
        GetShoreface(SS(n).MSLbeachWidths,SS(n).MHWbeachWidths,...
        SM(n).X1D,SM(n).Z1Dmean);
    fprintf('%i %6.3f %6.3f %6.4f\n',numel(ndx),SS(n).Slope,SlopeLF,rmseLF);
    ndx=find(SM(n).Z1Dmean >= 0.774 & SM(n).Z1Dmean <= 1.566);
    %fprintf('%i %6.3f %6.2f\n',numel(ndx),SS(n).Slope,nansum(SM(n).Z1Dmean(ndx)));
    %plot(-SM(n).X1D,SM(n).Z1Dmean)
    ss=[ss SlopeLF]; 
    plot(SS(n).Datenum,SlopeLF,'r+')
    hold on;plot(SS(n).Datenum,Slope,'k+')
    cdx=round(SS(n).MHWbeachWidths):round(SS(n).MSLbeachWidths);
    id1=find(round(SM(n).X1D) == round(SS(n).MHWbeachWidths)); 
    id2=find(round(SM(n).X1D) == round(SS(n).MSLbeachWidths)); 
    
    if ~isempty(cdx)
    plot([SS(n).Datenum SS(n).Datenum],[min(abs(gradient(SM(n).Z1Dmean(id1:id2))))...
        max(abs(gradient(SM(n).Z1Dmean(id1:id2))))],'m-')
    %plot(SS(n).Datenum,min(abs(gradient(SM(n).Z1Dmean(id1:id2)))),'mo')
    end
end
plot([SS(idx).Datenum],[SS(idx).Slope],'k-');
plot([SS(idx).Datenum],ss,'r-');

set(gca,'xlim',[-80,0],'ylim',[0.774 1.566])
