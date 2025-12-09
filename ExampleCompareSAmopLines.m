% Example to overlay a Mop range of Jumbo lines
close all
clearvars
Mop1=650%625%672;%625;%580;
Mop2=655%%674;%630;%590;

% load combined range of SA mop files
SA=SAcombineMops(Mop1,Mop2);


JumboIndexes=find(contains({SA.File}', 'jumbo','IgnoreCase',true)==1);

xmin=Inf;
figure('position',[22          54        1295         736]);
col=jet(5);
for n=1:5
ns=JumboIndexes(end-n+1);
i=find(SA(ns).Z > 0);
xmin=min([xmin SA(ns).X(i)']);
hold on;p(n)=plot(SA(ns).X(i),SA(ns).Y(i),'*','color',col(n,:),...
    'markersize',4,'linewidth',2,'DisplayName',datestr(SA(ns).Datenum));
end

for MopID=Mop1:Mop2
PlotLabelBackMopTransectUTM(MopID,'2d','k','ShadeOff')
end

set(gca,'dataaspectratio',[1 1 1],'fontsize',14)
xl=get(gca,'xlim');
set(gca,'xlim',[xmin-50 xl(2)],'color',[.7 .7 .7]);
legend(p,'location','westoutside','color',[.7 .7 .7])
xlabel('E UTM');ylabel('N UTM');
