load M00552SM.mat

figure('position',[204 99 1120 637]);
hold on;

% load BeachColorMap
% cm = BeachColorMap;
% 
idx=find([SM.Datenum] > datenum(2020,1,1) & strcmpi({SM.Source},'ig8wheel'));
for ns=idx%1:size(SM,2)
    
    n=length(SM(ns).X1D);
    x=SM(ns).X1D;
    y=SM(ns).Datenum*ones(1,n);
    %z=SM(ns).Z1Dmean;
    z=SM(ns).Z1Dtransect;
    p1=plot3(x,y,z,'-','linewidth',1);
    
    % if ~isempty(z)
    % x(isnan(z))=[];
    % y(isnan(z))=[];
    % z(isnan(z))=[]; 
    % x(z>4)=[];
    % y(z>4)=[];
    % z(z>4)=[]; 
    % end
    % 
    % if ~isempty(z)
    % y(x<0)=[];
    % z(x<0)=[];
    % x(x<0)=[];
    % %z(z>4)=NaN;  
    % 
    % zscaled = 1+size(BeachColorMap,1)*(z-BeachColorBarLims(1))/...
    % (BeachColorBarLims(2)-BeachColorBarLims(1));
    % zscaled(zscaled < 1)=1;
    % zscaled(zscaled > size(BeachColorMap,1))=size(BeachColorMap,1);   
    % 
    % scp=scatter3(x, y, z, 8, cm(ceil(zscaled),:), 'filled');
    % scp.MarkerFaceAlpha = .6;
    % scp.MarkerEdgeAlpha = .6;
    % 
    hold on;
    %end
end

view(3)
datetick(gca,'y');set(gca,'xdir','reverse')
grid on;

xl=get(gca,'xlim');yl=get(gca,'ylim');
fill3([xl xl(2) xl(1) xl(1)],[yl(1) yl yl(2) yl(1)],...
    [.774 .774 .774 .774 .774],[.8 .8 .8],'FaceAlpha',0.2);
xlabel('Xshore Distance (m)');
ylabel('Date');
zlabel('Elevation (m, NAVD88)');
title(['Mop ' num2str(SM(1).Mopnum) ' Transect Profiles']);
%title(['Mop ' num2str(SM(1).Mopnum) ' Area Mean Xshore Profiles']);
set(gca,'fontsize',13);
ytickangle(45);

% n=length(SM(1).X1D);
% x=SM(1).X1D;y=SM(1).Datenum*ones(1,n);z=SM(1).Z1Dmean;
% x(isnan(z))=[];y(isnan(z))=[];z(isnan(z))=[];   
% p1=plot3(x,y,z,'k-','linewidth',2);
% 
% n=length(SM(369).X1D);
% x=SM(369).X1D;y=SM(1).Datenum*ones(1,n);z=SM(369).Z1Dmean;
% x(isnan(z))=[];y(isnan(z))=[];z(isnan(z))=[];   
% p2=plot3(x,y,z,'k--','linewidth',2);
% legend([p1 p2],datestr(SM(1).Datenum),datestr(SM(369).Datenum),...
%     'location','northwest');

BeachColorbar

% pause
% 
% view(2)
% 
% pause
% 
% view(0,0)

mn=month(datetime([SM.Datenum],'convertfrom','datenum'));


figure('position',[224 79 1120 637]);
hold on;

idx=find(mn >= 11 | mn <= 5); 
for ns=idx  
    n=length(SM(ns).X1D);
    x=SM(ns).X1D;
    z=SM(ns).Z1Dtransect;
    if ~isempty(z)
    p1=plot(x,z,'r-','DisplayName','3 Winter Months (JFM)')   ; 
    end
    hold on;
end

idx=find(mn >= 6 & mn <= 11); 
for ns=idx  
    n=length(SM(ns).X1D);
    x=SM(ns).X1D;
    z=SM(ns).Z1Dtransect;
    if ~isempty(z)
    p2=plot(x,z,'b-','DisplayName','3 Summer Months (JAS)') ;  
    end
    hold on;
end

    % jumbo=find(contains({SM.File},'umbo') | contains({SM.File},'etski'));
    % 
    % x=SM(jumbo(end-2)).X1D;
    % z=SM(jumbo(end-2)).Z1Dtransect;
    % p3=plot(x,z,'c-','linewidth',3,'DisplayName',datestr(SM(jumbo(end-2)).Datenum)) ;
    % x=SM(jumbo(end-1)).X1D;
    % z=SM(jumbo(end-1)).Z1Dtransect;
    % p4=plot(x,z,'g-','linewidth',3,'DisplayName',datestr(SM(jumbo(end-1)).Datenum)) ; 
    % x=SM(jumbo(end)).X1D;
    % z=SM(jumbo(end)).Z1Dtransect;
    % p5=plot(x,z,'m-','linewidth',3,'DisplayName',datestr(SM(jumbo(end)).Datenum)) ; 
    
set(gca,'xdir','reverse')
grid on;

xl=get(gca,'xlim');yl=get(gca,'ylim');
% fill3([xl xl(2) xl(1) xl(1)],[yl(1) yl yl(2) yl(1)],...
%     [.774 .774 .774 .774 .774],[.8 .8 .8],'FaceAlpha',0.2);
xlabel('Xshore Distance (m)');
%ylabel('Date');
ylabel('Elevation (m, NAVD88)');
title(['Mop ' num2str(SM(1).Mopnum) ' Transect Profiles']);
%title(['Mop ' num2str(SM(1).Mopnum) ' Area Mean Xshore Profiles']);
set(gca,'fontsize',13);
legend([p1 p2 ],'location','northwest');
%legend([p1 p2 p3 p4 p5],'location','northwest');
%ytickangle(45);

