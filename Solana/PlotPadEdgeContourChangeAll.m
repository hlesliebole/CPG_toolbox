%  plot shorebox elevation changes all relative to a fixed xshore (Y) location
%  of a specified contour elevation (eg. the pad edge) on a specified 
%  starting survey date (16 apr 2024 for the solana pad).

% desired reference contour elevation for rectifying
zr=3.9; % approx contour of pad edge
%zr=2.5; % alternative of typical natural terrace edge contour

% xshore range to keep relative to the intial zr contour location
ymin=-100;
ymax=100;

% list of survey .mat files with gridded data in the solana shorebox
% coordinate frame which roughly translates to x=alongshore and y=xshore
% for this stretch of coast
sf=vertcat(...
'SolanaPostnourishmentGrid26Apr24.mat',...
'SolanaPostnourishmentGrid01May24.mat',...
'SolanaPostnourishmentGrid10May24.mat',...
'SolanaPostnourishmentGrid24May24.mat',...
'SolanaPostnourishmentGrid29May24.mat',...
'SolanaPostnourishmentGrid07Jun24.mat',...
'SolanaPostnourishmentGrid11Jun24.mat',...
'SolanaPostnourishmentGrid21Jun24.mat',...ls
'SolanaPostnourishmentGrid27Jun24.mat',...
'SolanaPostnourishmentGrid07Jul24.mat');

% load first postnourish shorebox grid X Y Z matrices
%   these all have dimensions ny x nx of the shorebox 1m resolution
%   grid
load SolanaPostnourishmentGrid16Apr24.mat % 16 apr 2024
% make initial baseline xshore rectified grid based on zr,ymin,ymax specs
%  Xr0,Yr0,Zr0 are the rectified and clipped matrices. zy = is the Y indices
%  of the zr contour location, which will be used to rectify all the other
%  surveys.
[Xr0,Yr0,Zr0,zy]=XshoreRectifyShoreboxBaseline(X,Y,Z,zr,ymin,ymax);

%% make a change figure
figure('position',[267    63   894   734]);
np=size(sf,1); % number of surveys to plot
nn=0;
% step through surveys. If first step, initialize with 16apr24 grid info
for n=1:np
    if n > 1
        Zr1=Zr2;
        iyzr1=iyzr2;
    else
        Zr1=Zr0;
        Xr1=Xr0;
        iyzr1=zy*0;
    end

% load next survey date X Y Z( matrices 
load(sf(n,:))
% rectify the survey relative to the initial 16 Apr zy location of the zr contour
%  iyzr2 are the Y indices of teh current location of the zr contour
[Xr2,Yr2,Zr2,iyzr2]=XshoreRectifyShorebox(X,Y,Z,zr,zy,ymin,ymax);

% make a left panel plot of the cumulative change
nn=nn+1;
subplot(np,2,nn)
% cumulative change of rectified grid is Zr2-Zr0
surf(Xr2,Yr2,Zr2-Zr0);colormap(flipud(polarmap));shading flat;%colorbar
view(2);set(gca,'clim',[-1.5 1.5],'ylim',[-20 20],'xlim',[775 2600])
hold on;
plot3(Xr1(1,:),iyzr2*0,iyzr2*0+10,'w-');
pl=plot3(Xr1(1,:),iyzr2,iyzr2*0+10,'k-');box on;grid on;
if nn < np*2-1;set(gca,'xticklabels',[]);
    if n == 1
        if zr == 3.9
         text(900,60,'Pad Edge (3.9m NAVD88 Contour) & Beach Face Change after 16 Apr PostNourishment Survey','fontsize',16);
         legend(pl,'Pad Edge','location','south')
        else
         text(500,60,'Natural Terrace Edge Elev (2.5m NAVD88 Contour) & Beach Face Change after 16 Apr PostNourishment Survey','fontsize',16);
         legend(pl,'Natural Terrac Elev','location','south')        
        end
        title('Cumulative Change','fontsize',16)
        
    end
    if n == 6
        if zr == 3.9
         ylabel('Xshore distance from Pad Edge Contour on 16 Apr (m)','fontsize',14)
        else
         ylabel('Xshore distance from Natural Terrace Edge Elev Contour on 16 Apr (m)','fontsize',14)
        end
    end
end
if n == np
        xlabel('S to N Alongshore Distance (m)','fontsize',14)
end

% make a right panel plot of the incremental change
nn=nn+1;
subplot(np,2,nn)
% Zr2-Zr1 is the increment rectified grid change
surf(Xr2,Yr2,Zr2-Zr1);colormap(flipud(polarmap));shading flat;%colorbar
view(2);set(gca,'clim',[-1.5 1.5],'ylim',[-20 20],'xlim',[775 2600])
hold on;
plot3(Xr1(1,:),iyzr2*0,iyzr2*0+10,'w-');
plot3(Xr1(1,:),iyzr2,iyzr1*0+10,'-','color',[0 0 0]);box on;grid on;
if nn < np*2-1;set(gca,'xticklabels',[]);end
ylabel([sf(n,26:30) ' '],'rotation',0,'fontsize',16)
if n == 1
        title('Incremental Change','fontsize',16)
end
if n == np
        xlabel('S to N Alongshore Distance (m)','fontsize',14)
end

end

% add colorbar at bottom
cb=colorbar('location','southoutside','fontsize',14);
cb.Position=[0.39    0.05501    0.2729    0.0054];
cb.Label.String='Elevation Change (m)';


%% supporting functions

function [Xr,Yr,Zr,zy]=XshoreRectifyShoreboxBaseline(X,Y,Z,zr,ymin,ymax)

% find xshore location of zr for each shorebox x

for nx=1:size(Z,2)
     iyz=find(Z(:,nx) < zr,1,'first');
    if ~isempty(iyz)
        zy(nx)=iyz;
    else
        zy(nx)=NaN;
    end
end

% make a Yr=ymin:ymax range xshore rectified version of the shorebox X Y Z 
%  where Yr=0 is the reference contour position 

Xr=NaN(numel(ymin:ymax),size(X,2));
Yr=NaN(numel(ymin:ymax),size(Y,2));
Zr=NaN(numel(ymin:ymax),size(Z,2));

for nx=1:size(Z,2)
    if ~isnan(zy(nx))
      Xr(:,nx)=X(zy(nx)+ymin:zy(nx)+ymax,nx);
      Zr(:,nx)=Z(zy(nx)+ymin:zy(nx)+ymax,nx);
      Yr(:,nx)=ymin:ymax;
    end
end

end

%% 

function [Xr,Yr,Zr,iyzr]=XshoreRectifyShorebox(X,Y,Z,zr,zy,ymin,ymax)

% find xshore location of zr for each shorebox x

for nx=1:size(Z,2)
     iyz=find(Z(:,nx) < zr,1,'first');
    if ~isempty(iyz)
        iyzr(nx)=iyz-zy(nx);
    else
        iyzr(nx)=NaN;
    end
end

% make a Yr=ymin:ymax range xshore rectified version of the shorebox X Y Z 
%  where Yr=0 is the reference contour position 

Xr=NaN(numel(ymin:ymax),size(X,2));
Yr=NaN(numel(ymin:ymax),size(Y,2));
Zr=NaN(numel(ymin:ymax),size(Z,2));

for nx=1:size(Z,2)
    if ~isnan(zy(nx))
      Xr(:,nx)=X(zy(nx)+ymin:zy(nx)+ymax,nx);
      Zr(:,nx)=Z(zy(nx)+ymin:zy(nx)+ymax,nx);
      Yr(:,nx)=ymin:ymax;
    end
end

end