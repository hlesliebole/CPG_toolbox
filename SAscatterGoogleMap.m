load MopTableUTM.mat

MopStart=find(strcmp([Mop.Name],'D0645'));
MopEnd=find(strcmp([Mop.Name],'D0657'));
MopStart=find(strcmp([Mop.Name],'D0632'));
MopEnd=find(strcmp([Mop.Name],'D0647'));
MopStart=find(strcmp([Mop.Name],'D0644'));
MopEnd=find(strcmp([Mop.Name],'D0665'));

CS=SAcombineMops(MopStart,MopEnd);
idx=find(strcmp({CS.Source},'Multibeam'));
%idx=find([CS.Datenum] == datenum(2023,11,21));
%idx=idx(2);
idx=idx(end);

x=vertcat(CS(idx).X);
y=vertcat(CS(idx).Y);
z=vertcat(CS(idx).Z);
% trim higher elevations
zidx=find(z < 6);
x=x(zidx);y=y(zidx);z=z(zidx);
% convert to lat lon
[lat,lon]=utm2deg(x,y,repmat('11 S',[length(x) 1]));

%figure;
figure('position',[1          12        1275         783]);%plot(lon,lat,'y.');
for n=MopStart:MopEnd
plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],'m-')
text(Mop.OffLon(n),Mop.OffLat(n),num2str(n),'horizontalalign','right','color','w')
end
hold on;
[ScatterPlot,ColorBarPlot]=ColorScatter2d(lon,lat,z);
hold on
load MopTableUTM.mat
for n=MopStart:MopEnd
plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],'m-')
text(Mop.OffLon(n),Mop.OffLat(n),num2str(n),'horizontalalign','right','color','w')
end
%plot(lon2(idx2),lat2(idx2),'m.')
plot_google_map('MapType', 'satellite')
set(gca,'clim',[-13 -4])
set(gca,'fontsize',16)

title({[datestr(CS(idx).Datenum) ' |  Multibeam '],...
    'Mops 645 to 655'},...
    'fontsize',16)

makepng(['Multibeam' datestr(CS(idx).Datenum,'YYYYmmDD') 'Mops' num2str(MopStart) 'to' num2str(MopEnd) '.png'])