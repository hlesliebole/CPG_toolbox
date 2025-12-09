% plots the time series of the mode of a Mop areas lidar elevation

MopNumber=625;
MinModeElev=2; % look above 2m navd88
Zres=0.1; % elevation bin size prior to mode calc

% load the SA.mat file
SAmatfile=['M' num2str(MopNumber,'%5.5i') 'SA.mat'];
fprintf('Loading SA struct array from %s\n',SAmatfile)
load(SAmatfile);

fprintf('\nMop %i LiDAR Data Summary\n\n',MopNumber);

% % truck, atv and airborne lidar
% ldx=find( ( strcmp({SA.Source},'Trk') | strcmp({SA.Source},'AtvMR') | ...
%     strcmp({SA.Source},'TrkMR') | strcmp({SA.Source},'UTAir') | strcmp({SA.Source},'KMair')) );

% just truck and atv
ldx=find( ( strcmp({SA.Source},'Trk') | strcmp({SA.Source},'AtvMR') | ...
    strcmp({SA.Source},'TrkMR') ) );

fprintf(' 1. %i total LiDAR surveys\n',numel(ldx));
fprintf(' 2. First Survey %s \n',datetime(SA(ldx(1)).Datenum,'convertfrom','datenum'));
fprintf(' 3. Last Survey %s \n',datetime(SA(ldx(end)).Datenum,'convertfrom','datenum'));

% reduce SA to just the lidar surveys
SA=SA(ldx);

% loop through surveys and plot mode of each above the min mode elevation

Sdatetime=datetime([],[],[]);
Zmode=[];
for n=1:size(SA,2)
    z=Zres*round([SA(n).Z]/Zres);
    z(z < MinModeElev) =[];
    Sdatetime(n)=datetime(SA(n).Datenum,'convertfrom','datenum');
    Zmode(n)=mode(z);
end

figure('position',[100,200,900,400]);
plot(Sdatetime,Zmode,'*');grid on;set(gca,'fontsize',16)
ylabel('Mode(Z) (m, navd88)')
title('Mode of Lidar Elevations above 2m navd88 (0.1m Vertical Resolution)');
