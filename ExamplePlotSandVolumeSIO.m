% Example code to make a plot of the 1m gridded data from a Jumbo
%  survey for a combined range of mop SG mat file in a new CG struct array.

%% Uses the m-script function CG=CombineSGdata(mapth,MopStart,MopEnd)

%% settings
%mpath='/volumes/group/MOPS/'; % path to cpg mop files
mpath='/Users/William/Desktop/MOPS/';
MopStart=507%666;%582; % start mop number
MopEnd=514%683;%590; % ending mop number

%% load the combined SG mat file data.  Return the combined
%  data to a struct array SG instead of CG to use the normal 
%  single mop SG gridding and plotting code with the combined data.
SG=CombineSGdata(mpath,MopStart,MopEnd);

%% identify jumbos by checking the original survey file names for the 
%  word jumbo
utair=find(contains({SG.Source},'UTAir'));
trk=find(contains({SG.Source},'Trk'));
usace=find(contains({SG.Source},'USACE'));
usgs=find(contains({SG.Source},'USGS'));
ccc=find(contains({SG.Source},'CCC'));
utrk=unique([utair trk usace usgs ccc]);

% datestr([SG(jumbo).Datenum])
% fprintf('The SG struct array has %i UTAir Surveys.\n',numel(jumbo))
% 
% fprintf('Plotting the most recent airborne Survey in the database.\n',numel(jumbo))

%% Need to turn the saved SG struct array grid points into a 2d grid array

% make grid area encompassing all survey data
minx=min(vertcat(SG.X));
maxx=max(vertcat(SG.X));
miny=min(vertcat(SG.Y));
maxy=max(vertcat(SG.Y));
%[X,Y]=meshgrid(min(x):max(x),min(y):max(y));
[X,Y]=meshgrid(minx:maxx,miny:maxy);

% inlet location
[e,n,zn]=deg2utm(32.9338,-117.26);

figure('position',[ 21          56        1385         728]);

%SurvNum=trk(end-1);
SurvNum=size(SG,2);
Zmn=X*0+100;
for SurvNum=1:size(SG,2)
% Mop area 1m grid points with valid data
x=SG(SurvNum).X;
y=SG(SurvNum).Y;
% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);
Z1=X*NaN; % initialize the 2d elevation Z array as NaNs
Z1(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices
Zmn=min(Zmn,Z1);
end
Zmn(Zmn == 100)=NaN;

% now find subaerial volumes
mm=0;
v=[];
nv=[];
for SurvNum=utrk(2:end)
    mm=mm+1;
% Mop area 1m grid points with valid data
x=SG(SurvNum).X;
y=SG(SurvNum).Y;
% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);
Z1=X*NaN; % initialize the 2d elevation Z array as NaNs
Z1(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices
Z1(Z1 < 0.774)=NaN;%Z1(Z1 < 1.344)=NaN;%
Z1(Z1 > 3.0)=NaN;

dz=Z1-Zmn;
v(mm)=sum(dz(:),'omitnan');
if v(mm) < 0
    SurvNum
end
nv(mm)=numel(find(~isnan(dz(:))));
end

v(nv < 100)=NaN;

plot([SG(utrk(2:end)).Datenum],v/1000,'k.','markersize',35);datetick
grid on;
set(gca,'fontsize',18);
xlabel('Date');
ylabel('Mobile Sand Volume (m^{3} x 1000)');
title({'Measured Volume of "Known to be Mobile" Beach Sand',...
    '(Sand Shoreward of MSL) on SIO''s 400m Beachfront'},...
    'fontsize',24);
box on
set(gca,'linewidth',2)
makepng('SIOsandVolume.png')


%figure
%Z1=dz;
% plot the xutm,yutm,z 3d surface

% %ax1=axes('position',[0.05    0.0500    0.1504    0.95]);
% p1=surf(X,Y,Z1,'linestyle','none');
% view(2)
% %BeachColorbar
% axis equal
% hold on;plot(e,n,'k.','markersize',20);
% 
% %set(gca,'xlim',[473300 474100]);hold on;
% set(gca,'dataaspectratio',[1 1 1])
% 
% grid on;
% set(gca,'color',[.7 .7 .7],'fontsize',14);
% 
% y_labels = get(gca, 'YTick');
% set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
% ylabel('northings (m)');
% x_labels = get(gca, 'XTick');
% set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
% xlabel('eastings (m)');
% zlabel('Elevation (m.NAVD88)');
% 
% view(2)
% %axis equal
% 
% d=Z2-Z1;erode=round(sum(d(d<0))/1000);dep=round(sum(d(d>0))/1000);
% str=['-V= ' num2str(erode) 'K  ; +V= ' num2str(dep) 'K'];str
% text(473350,3653120,str,'fontweight','bold','fontsize',18);
% 
% set(gca,'clim',[-3 3]);
% polarmap;
% 
% title({
%     datestr(SG(SurvNum).Datenum), ' to ', datestr(SG(SurvNum2).Datenum)},...
%     'fontsize',16)
% 
% %% -------------------------------------
% %% 2017 
% 
% % first survey
% 
% %didx=find([SG(jumbo).Datenum] ==  datenum(2016,9,29)); % Cardiff
% didx=find([SG(jumbo).Datenum] ==  datenum(2016,9,27)); % Torrey
% %didx=find([SG(jumbo).Datenum] ==  datenum(2016,7,21)); % Torrey
% 
% 
% SurvNum=jumbo(didx); 
% % Mop area 1m grid points with valid data
% x=SG(SurvNum).X;
% y=SG(SurvNum).Y;
% 
% % Make 2d x,y utm arrays encompassing the valid points
% idx=sub2ind(size(X),y-miny+1,x-minx+1);
% 
% Z1=X*NaN; % initialize the 2d elevation Z array as NaNs
% Z1(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices
% 
% % nextjumbo survey
% % didx=find([SG(jumbo).Datenum] ==  datenum(2017,1,26)); % Cardiff
% didx=find([SG(jumbo).Datenum] ==  datenum(2017,4,11)); % Torrey
% 
% SurvNum2=jumbo(didx);
% 
% x=SG(SurvNum2).X;
% y=SG(SurvNum2).Y;
% 
% % Make 2d x,y utm arrays encompassing the valid points
% idx=sub2ind(size(X),y-miny+1,x-minx+1);
% 
% Z2=X*NaN; % initialize the 2d elevation Z array as NaNs
% Z2(idx)=SG(SurvNum2).Z; % overlay the valid data points using the 1d indices
% 
% % plot the xutm,yutm,z 3d surface
% 
% ax2=axes('position',[0.35    0.0500    0.1504    0.95]);
% p2=surf(X,Y,Z2-Z1,'linestyle','none');
% hold on;plot(e,n,'k.','markersize',20);
% 
% %set(gca,'xlim',[473300 474100]);hold on;
% set(gca,'dataaspectratio',[1 1 1])
% 
% grid on;
% set(gca,'color',[.7 .7 .7],'fontsize',14);
% 
% y_labels = get(gca, 'YTick');
% set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
% ylabel('northings (m)');
% x_labels = get(gca, 'XTick');
% set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
% xlabel('eastings (m)');
% zlabel('Elevation (m.NAVD88)');
% 
% view(2)
% %axis equal
% 
% d=Z2-Z1;erode=round(sum(d(d<0))/1000);dep=round(sum(d(d>0))/1000);
% str=['-V= ' num2str(erode) 'K  ; +V= ' num2str(dep) 'K'];str
% text(473350,3653120,str,'fontweight','bold','fontsize',18);
% 
% set(gca,'clim',[-3 3]);
% polarmap;
% 
% title({
%     datestr(SG(SurvNum).Datenum), ' to ', datestr(SG(SurvNum2).Datenum)},...
%     'fontsize',16)
% 
% %% -------------------------------------
% 
% %% -------------------------------------
% %% 2021 
% 
% % first survey
% 
% %didx=find([SG(jumbo).Datenum] ==  datenum(2020,10,1)); % Cardiff
% didx=find([SG(jumbo).Datenum] ==  datenum(2020,10,14)); % Torrey
% 
% SurvNum=jumbo(didx); 
% % Mop area 1m grid points with valid data
% x=SG(SurvNum).X;
% y=SG(SurvNum).Y;
% 
% % Make 2d x,y utm arrays encompassing the valid points
% idx=sub2ind(size(X),y-miny+1,x-minx+1);
% 
% Z1=X*NaN; % initialize the 2d elevation Z array as NaNs
% Z1(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices
% 
% % nextjumbo survey
% %didx=find([SG(jumbo).Datenum] ==  datenum(2021,1,28)); % Cardiff
% didx=find([SG(jumbo).Datenum] ==  datenum(2021,2,10)); % Torrey
% SurvNum2=jumbo(didx);
% 
% x=SG(SurvNum2).X;
% y=SG(SurvNum2).Y;
% 
% % Make 2d x,y utm arrays encompassing the valid points
% idx=sub2ind(size(X),y-miny+1,x-minx+1);
% 
% Z2=X*NaN; % initialize the 2d elevation Z array as NaNs
% Z2(idx)=SG(SurvNum2).Z; % overlay the valid data points using the 1d indices
% 
% % plot the xutm,yutm,z 3d surface
% 
% ax3=axes('position',[0.35    0.0500    0.1504    0.95]);
% p2=surf(X,Y,Z2-Z1,'linestyle','none');
% hold on;plot(e,n,'k.','markersize',20)
% 
% %set(gca,'xlim',[473300 474100]);hold on;
% set(gca,'dataaspectratio',[1 1 1])
% 
% grid on;
% set(gca,'color',[.7 .7 .7],'fontsize',14);
% 
% y_labels = get(gca, 'YTick');
% set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
% ylabel('northings (m)');
% x_labels = get(gca, 'XTick');
% set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
% xlabel('eastings (m)');
% zlabel('Elevation (m.NAVD88)');
% 
% view(2)
% %axis equal
% 
% d=Z2-Z1;erode=round(sum(d(d<0))/1000);dep=round(sum(d(d>0))/1000);
% str=['-V= ' num2str(erode) 'K  ; +V= ' num2str(dep) 'K'];str
% text(473350,3653120,str,'fontweight','bold','fontsize',18);
% 
% 
% 
% set(gca,'clim',[-3 3]);
% polarmap;
% 
% title({
%     datestr(SG(SurvNum).Datenum), ' to ', datestr(SG(SurvNum2).Datenum)},...
%     'fontsize',16)
% 
% %% -------------------------------------
% %% -------------------------------------
% %% 2023 
% 
% % first survey
% 
% %didx=find([SG(jumbo).Datenum] ==  datenum(2022,10,12)); % Cardiff
% didx=find([SG(jumbo).Datenum] ==  datenum(2022,10,10)); % Torrey
% 
% 
% SurvNum=jumbo(didx); 
% % Mop area 1m grid points with valid data
% x=SG(SurvNum).X;
% y=SG(SurvNum).Y;
% 
% % Make 2d x,y utm arrays encompassing the valid points
% idx=sub2ind(size(X),y-miny+1,x-minx+1);
% 
% Z1=X*NaN; % initialize the 2d elevation Z array as NaNs
% Z1(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices
% 
% % nextjumbo survey
% %didx=find([SG(jumbo).Datenum] ==  datenum(2023,1,23)); % Cardiff
% didx=find([SG(jumbo).Datenum] ==  datenum(2023,1,24)); % Torrey
% SurvNum2=jumbo(didx);
% 
% x=SG(SurvNum2).X;
% y=SG(SurvNum2).Y;
% 
% % Make 2d x,y utm arrays encompassing the valid points
% idx=sub2ind(size(X),y-miny+1,x-minx+1);
% 
% Z2=X*NaN; % initialize the 2d elevation Z array as NaNs
% Z2(idx)=SG(SurvNum2).Z; % overlay the valid data points using the 1d indices
% 
% % plot the xutm,yutm,z 3d surface
% 
% ax4=axes('position',[0.35    0.0500    0.1504    0.95]);
% p2=surf(X,Y,Z2-Z1,'linestyle','none');
% hold on;plot(e,n,'k.','markersize',20)
% 
% %set(gca,'xlim',[473300 474100]);hold on;
% set(gca,'dataaspectratio',[1 1 1])
% 
% grid on;
% set(gca,'color',[.7 .7 .7],'fontsize',14);
% 
% y_labels = get(gca, 'YTick');
% set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
% ylabel('northings (m)');
% x_labels = get(gca, 'XTick');
% set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
% xlabel('eastings (m)');
% zlabel('Elevation (m.NAVD88)');
% 
% view(2)
% %axis equal
% 
% d=Z2-Z1;erode=round(sum(d(d<0))/1000);dep=round(sum(d(d>0))/1000);
% str=['-V= ' num2str(erode) 'K  ; +V= ' num2str(dep) 'K'];str
% text(473350,3653120,str,'fontweight','bold','fontsize',18);
% 
% set(gca,'clim',[-3 3]);
% polarmap;
% 
% % title({
% %     datestr(SG(SurvNum).Datenum), ' to ', datestr(SG(SurvNum2).Datenum)},...
% %     'fontsize',16)
% 
% title({
%     datestr(SG(SurvNum).Datenum), ' to ','23-Jan-2023'},...
%     'fontsize',16)
% 
% for m=578:595
% plot3([Mop.BackXutm(m) Mop.OffXutm(m)],[Mop.BackYutm(m) Mop.OffYutm(m)],[10 10],'m.-',...
%     'markersize',15,'linewidth',1);
% text(Mop.BackXutm(m)+50,Mop.BackYutm(m),num2str(m),'color','w','fontsize',14,...
%     'backgroundcolor','k')
% end
% 
% %t=text(ax4,473550,3652000,10,{'PLEASE','STAND BY'},'fontweight','bold','fontsize',20)
% %% -------------------------------------
% 
% cb=colorbar;
% cb.Label.String='Elevation Change (m)';
% cb.FontSize=14;
% 
% set(ax1,'position',[0.075    0.0500    0.1504    0.95]);
% set(ax2,'position',[0.3    0.0500    0.1504    0.95]);
% set(ax3,'position',[0.525    0.0500    0.1504    0.95]);
% set(ax4,'position',[0.75    0.0500    0.1504    0.95]);
% colormap(flipud(colormap))
% 
% makepng('TorreyNorthSevereWinterChanges.png')