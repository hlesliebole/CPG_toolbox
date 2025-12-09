ztide=[-0.631 -.058 .218 .774 1.344 1.566 2.119 ];

figure
idx=kmeans(slope(:),6);ScatterPlotMopUTM(X(:),Y(:),idx,'2d');

%slp=slope(:);
figure
for n=min(idx):max(idx)
    plot(1,min(slope(idx == n)),'k*')
    hold on
    plot([1 2],[min(slope(idx == n)) min(zg(idx == n))],'k:')
    hold on
end

for n=1:length(ztide)
    plot(2,ztide(n),'k+')    
end

set(gca,'xlim',[0 3]);grid on;
text(2.1,ztide(1),'LAT');
text(2.1,ztide(2),'MLLW');
text(2.1,ztide(3),'MLW');
text(2.1,ztide(4),'MSL');
text(2.1,ztide(5),'MHW');
text(2.1,ztide(6),'MHHW');
text(2.1,ztide(7),'HAT');



