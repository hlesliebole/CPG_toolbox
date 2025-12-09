function [fh]=PlotCMstruct(varargin)
 
% up to two inputs: 
%     PlotCMstruct([CMstructName])
%     PlotCMstruct([CMstructName],[SurvNum])

%  Plots the data is the CPGMOP structural array MopStruct. 
%  If a SurveyNum is also input, it just plots the data for
%  that struct array index number.
%
%  returns the figure handle 

% assign variables

% figure out if it is a single MOP struct array or a nested set of
%  MOP struct arrays
S=varargin{1};
if size(fieldnames(varargin{1}),1) > 1  % single MOP struct array 
MopStructType=inputname(1); % get the actual input struct variable name
isnest=0; % no nest logical
nmops=1;
MopNumber=S(1).Mopnum;
else % already a nested struct of mops
MopStructType=fieldnames(S); % Mop struct type is nested field name
isnest=1; % yes nest logical
nmops=size(S,2);
MopNumber=S(1).SA(1).Mopnum;
end

% assign the survey datenum to plot if one was input
if nargin == 2;SurvNum=varargin{2};else;SurvNum=[];end    


%-------------------------------------------------------------------
%   SA  survey data
%-------------------------------------------------------------------

if strcmp(MopStructType,'SA') 
    
    for m=1:nmops % loop through Mop areas
       
         if isempty(SurvNum)
            if isnest
             x=vertcat(S(m).SA(1).X);
             y=vertcat(S(m).SA(1).Y);
             z=vertcat(S(m).SA(1).Z);
            else
             x=vertcat(S(1).X);
             y=vertcat(S(1).Y);
             z=vertcat(S(1).Z);
            end
         else
            if isnest
             x=vertcat(S(m).SA(SurvNum).X);
             y=vertcat(S(m).SA(SurvNum).Y);
             z=vertcat(S(m).SA(SurvNum).Z);
            else
             tidx=find(vertcat(S.Datenum) == SurvNum);
             x=vertcat(S(SurvNum).X);
             y=vertcat(S(SurvNum).Y);
             z=vertcat(S(SurvNum).Z);
            end
         end
           
            % screen for NaN elevations
            x(isnan(z))=[];y(isnan(z))=[];z(isnan(z))=[];
            % screen for extreme min-max elevations
            x(z < -15)=[];y(z < -15)=[];z(z < -15)=[];
            x(z > 10)=[];y(z > 10)=[];z(z > 10)=[];
            PlotMethod='3d';
            %PlotMopTransectUTM(MopNumber-2,PlotMethod,'k','ShadeOff');
            PlotMopTransectUTM(MopNumber-1,PlotMethod,'k','ShadeOff');
            PlotMopTransectUTM(MopNumber,PlotMethod,'m','ShadeOn');
            PlotMopTransectUTM(MopNumber+1,PlotMethod,'k','ShadeOff');
            %PlotMopTransectUTM(MopNumber+2,PlotMethod,'k','ShadeOff');
     
            ScatterPlotBeachUTM(x,y,z,'3Dcolor');
            set(gca,'color',[.6 .6 .6]);
             hold on;
             view(3)

    end
    
 grid on;
 BeachColorbar 
end

%-------------------------------------------------------------------
%   SG  gridded  data
%-------------------------------------------------------------------

if strcmp(MopStructType,'SG') 
    
            PlotMethod='3d';
            %PlotMopTransectUTM(MopNumber-2,PlotMethod,'k','ShadeOff');
            PlotMopTransectUTM(MopNumber-1,PlotMethod,'k','ShadeOff');
            PlotMopTransectUTM(MopNumber,PlotMethod,'m','ShadeOn');
            PlotMopTransectUTM(MopNumber+1,PlotMethod,'k','ShadeOff');
            %PlotMopTransectUTM(MopNumber+2,PlotMethod,'k','ShadeOff');
     
            
x=S(SurvNum).X;
y=S(SurvNum).Y;


[X,Y]=meshgrid(min(x):max(x),min(y):max(y));
idx=sub2ind(size(X),y-min(y)+1,x-min(x)+1);
Z=X*NaN;

% plot surface
Z(idx)=S(SurvNum).Z;
p1=surf(X,Y,Z,'linestyle','none');
    
grid on;
set(gca,'color',[.7 .7 .7],'fontsize',14);


y_labels = get(gca, 'YTick');
set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
ylabel('northings (m)');
x_labels = get(gca, 'XTick');
set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
xlabel('eastings (m)');


 BeachColorbar 
end 

%-------------------------------------------------------------------
%   CG  gridded  data
%-------------------------------------------------------------------

if strcmp(MopStructType,'CG') 
    
            PlotMethod='3d';
            for MopNumber=min([S(1).Mopnum]):max([S(1).Mopnum])
            %PlotMopTransectUTM(MopNumber-2,PlotMethod,'k','ShadeOff');
            %PlotMopTransectUTM(MopNumber-1,PlotMethod,'k','ShadeOff');
            %PlotMopTransectUTM(MopNumber,PlotMethod,'m','ShadeOn');
            PlotLabelMopTransectUTM(MopNumber,PlotMethod,'k','ShadeOff');
            %PlotMopTransectUTM(MopNumber+2,PlotMethod,'k','ShadeOff');
            %end
            end
     
            
x=S(SurvNum).X;
y=S(SurvNum).Y;


[X,Y]=meshgrid(min(x):max(x),min(y):max(y));
idx=sub2ind(size(X),y-min(y)+1,x-min(x)+1);
Z=X*NaN;

% plot surface
Z(idx)=S(SurvNum).Z;
p1=surf(X,Y,Z,'linestyle','none');
    
grid on;
set(gca,'color',[.7 .7 .7],'fontsize',14);


y_labels = get(gca, 'YTick');
set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
ylabel('northings (m)');
x_labels = get(gca, 'XTick');
set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
xlabel('eastings (m)');


BeachColorbar; 
end 


%-------------------------------------------------------------------
%   SM  morpho data
%-------------------------------------------------------------------

if strcmp(MopStructType,'SM') 
    
    for m=1:nmops % loop through Mop areas
         if isempty(SurvNum)
            x=vertcat(S(m).SM(:).X1D);
            zt=vertcat(S(m).SM(:).Z1Dtransect);
            zmn=vertcat(S(m).SM(:).Z1Dmean);
            zmd=vertcat(S(m).SM(:).Z1Dmedian);
            zstd=vertcat(S(m).SM(:).Z1Dstd);
            zmin=vertcat(S(m).SM(:).Z1Dmin);
            zmax=vertcat(S(m).SM(:).Z1Dmax);
         else
            tidx=SurvNum;
            x=vertcat(S(tidx).X1D);
            zt=vertcat(S(tidx).Z1Dtransect);
            zmn=vertcat(S(tidx).Z1Dmean);
            zmd=vertcat(S(tidx).Z1Dmedian);
            zstd=vertcat(S(tidx).Z1Dstd);
            zmin=vertcat(S(tidx).Z1Dmin);
            zmax=vertcat(S(tidx).Z1Dmax);
         end
            
            % reduce to profile section with valid data
            igood=find(~isnan(zmn));
            x=x(igood);
            zt=zt(igood);
            zmn=zmn(igood);
            zmd=zmd(igood);
            zstd=zstd(igood);
            zmin=zmin(igood);
            zmax=zmax(igood);
            
            set(gca,'linewidth',2);hold on;
            p1=plot(x,zt,'k-','linewidth',2);set(gca,'xdir','reverse','fontsize',14);
            hold on;
            p2=plot(x,zmn,'g-','linewidth',2);set(gca,'xdir','reverse','fontsize',14);
            p3=plot(x,zmn,'b:','linewidth',2);set(gca,'xdir','reverse','fontsize',14);
            p4=plot(x,zmin,'r-','linewidth',1);set(gca,'xdir','reverse','fontsize',14);
            p5=plot(x,zmax,'r-','linewidth',1);set(gca,'xdir','reverse','fontsize',14);
            
            % label tide elevations
            ztide=[-0.631 -.058 .218 .774 1.344 1.566 2.119 ];
            set(gca,'xlim',[x(find(~isnan(zmn), 1 )) x(find(~isnan(zmn), 1, 'last' ))]); 
            xl=get(gca,'xlim');
            for n=1:length(ztide)
                if n == 4
                 plot(xl,[ztide(n) ztide(n)],'b--');  
                else
                 plot(xl,[ztide(n) ztide(n)],'b:','linewidth',1);  
                end
            end

        text(xl(1),ztide(1),' LAT','fontsize',14);
        text(xl(1),ztide(2),' MLLW');
        text(xl(1),ztide(3),' MLW');
        text(xl(1),ztide(4),' MSL','fontsize',14);
        text(xl(1),ztide(5),' MHW');
        text(xl(1),ztide(6),' MHHW');
        text(xl(1),ztide(7),' HAT','fontsize',14);
            
        legend([p1 p2 p3 p4],'Transect Profile (from gridded data)','Area Mean Profile',...
                'Area Median Profile','Area Min-Max Profile Elevations','location','northwest');
            
            
        ylabel('Elevation (m, NAVD88)');
        xlabel('Xshore Distance from MOP Back Beach Point (m)');
        hold on;

    end

 grid on;
 
%  % fix axes based on global morpho limits for all surveys
%  set(gca,'xlim',[min([S.X1D]) max([S.X1D])],...
%      'ylim',[min([S.Z1Dmean]) max([S.Z1Dmean])]);

end 

%------------------------------- 
% Global Morphology 1D Profiles
%------------------------------- 
if strcmp(MopStructType,'GM') && SurvNum == 1 
                    
    set(gca,'linewidth',2);hold on;
    p1=plot(S.X1D,S.Z1Dtransect,'k-','linewidth',2);set(gca,'xdir','reverse','fontsize',14);
    hold on;
    p2=plot(S.X1D,S.Z1Dmean,'g-','linewidth',2);set(gca,'xdir','reverse','fontsize',14);
    p3=plot(S.X1D,S.Z1Dmedian,'g:','linewidth',2);set(gca,'xdir','reverse','fontsize',14);
    p4=plot(S.X1D,S.Z1Dmin,'r-','linewidth',1);set(gca,'xdir','reverse','fontsize',14);
    p5=plot(S.X1D,S.Z1Dmax,'r-','linewidth',1);set(gca,'xdir','reverse','fontsize',14);
            
    ztide=[-0.631 -.058 .218 .774 1.344 1.566 2.119 ];
    xl=get(gca,'xlim');
     for n=1:length(ztide)
       if n == 4
           plot(xl,[ztide(n) ztide(n)],'b--');  
       else
           plot(xl,[ztide(n) ztide(n)],'b:','linewidth',1);  
       end
     end

     text(xl(1),ztide(1),' LAT','fontsize',14);
     text(xl(1),ztide(2),' MLLW');
     text(xl(1),ztide(3),' MLW');
     text(xl(1),ztide(4),' MSL','fontsize',14);
     text(xl(1),ztide(5),' MHW');
     text(xl(1),ztide(6),' MHHW');
     text(xl(1),ztide(7),' HAT','fontsize',14);
            
     legend([p1 p2 p3 p4],'Global Mean Transect Profile','Global Area Mean Profile',...
            'Global Area Median Profile',...
                'Global Area Min-Max Profile Elevations','location','northwest');
                   
     ylabel('Elevation (m, NAVD88)');
     xlabel('Xshore Distance from MOP Back Beach Point (m)');
     hold on;

     grid on;

end 

%------------------------------- 
% Global Morphology 2D Bathy
%------------------------------- 
if strcmp(MopStructType,'GM')  && SurvNum == 2 

           
            PlotMethod='3d';
            %PlotMopTransectUTM(MopNumber-2,PlotMethod,'k','ShadeOff');
            PlotMopTransectUTM(MopNumber-1,PlotMethod,'k','ShadeOff');
            PlotMopTransectUTM(MopNumber,PlotMethod,'m','ShadeOn');
            PlotMopTransectUTM(MopNumber+1,PlotMethod,'k','ShadeOff');
            %PlotMopTransectUTM(MopNumber+2,PlotMethod,'k','ShadeOff');
               
            % convert gridded x,y,z points to 2d matrices
x=S.X2D;
y=S.Y2D;


[X,Y]=meshgrid(min(x):max(x),min(y):max(y));
idx=sub2ind(size(X),y-min(y)+1,x-min(x)+1);
Z=X*NaN;

% plot surfaces 
Z(idx)=S.Z2Dmin;
p1=surf(X,Y,Z,'linestyle','none','facecolor',[.5 .5 .5],'facealpha',.9);
hold on;
Z(idx)=S.Z2Dmean;
p2=surf(X,Y,Z,'linestyle','none','facecolor','b','facealpha',.5);
Z(idx)=S.Z2Dmedian;
p3=surf(X,Y,Z,'linestyle','none','facecolor','g','facealpha',.5);
Z(idx)=S.Z2Dmax;
p4=surf(X,Y,Z,'linestyle','none','facecolor','r','facealpha',.5);
legend([p1 p2 p3 p4],'Global Min','Global Mean','Global Median',...
    'Global Max','location','northwest')
    
grid on;
set(gca,'color',[.7 .7 .7],'fontsize',14);
 
y_labels = get(gca, 'YTick');
set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
ylabel('Northings (m)');
x_labels = get(gca, 'XTick');
set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
xlabel('Eastings (m)');
 
 %BeachColorbar 
end 

%------------------------------- 
% Global Monthly Morphology 2D Bathy
%------------------------------- 

if strcmp(MopStructType,'GM')  && SurvNum > 2
    
    mon=SurvNum-2;

           
            PlotMethod='3d';
            %PlotMopTransectUTM(MopNumber-2,PlotMethod,'k','ShadeOff');
            PlotMopTransectUTM(MopNumber-1,PlotMethod,'k','ShadeOff');
            PlotMopTransectUTM(MopNumber,PlotMethod,'m','ShadeOn');
            PlotMopTransectUTM(MopNumber+1,PlotMethod,'k','ShadeOff');
            %PlotMopTransectUTM(MopNumber+2,PlotMethod,'k','ShadeOff');
               
            % convert gridded x,y,z points to 2d matrices
x=S.X2D;
y=S.Y2D;

[X,Y]=meshgrid(min(x):max(x),min(y):max(y));
idx=sub2ind(size(X),y-min(y)+1,x-min(x)+1);
Z=X*NaN;

%plot global min surface 
Z(idx)=S.Z2Dmin;
p1=surf(X,Y,Z,'linestyle','none','facecolor',[.5 .5 .5],'facealpha',.9);

% month surface
x=S.MM(mon).X2D;
y=S.MM(mon).Y2D;

[X,Y]=meshgrid(min(x):max(x),min(y):max(y));
idx=sub2ind(size(X),y-min(y)+1,x-min(x)+1);
Z=X*NaN;

% plot surfaces 

hold on;
Z(idx)=S.MM(mon).Z2Dmean;
p2=surf(X,Y,Z,'linestyle','none','facecolor','b','facealpha',.5);
Z(idx)=S.MM(mon).Z2Dmedian;
p3=surf(X,Y,Z,'linestyle','none','facecolor','g','facealpha',.5);

legend([p1 p2 p3],'Global Min','Global Month Mean','Global Month Median',...
    'location','northwest')
    
grid on;
set(gca,'color',[.7 .7 .7],'fontsize',14);
 
y_labels = get(gca, 'YTick');
set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
ylabel('Northings (m)');
x_labels = get(gca, 'XTick');
set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
xlabel('Eastings (m)');
 

end 

%------------------------------- 
% Wave Parameters for all surveys
%------------------------------- 

if strcmp(MopStructType,'SW') 

    subplot(3,1,1)
    set(gca,'fontsize',14);hold on;
%     p1=plot([S.Datenum],[S.Emean],'ks');
%     p2=plot([S.Datenum],[S.Emedian],'bo');
%     p3=plot([S.Datenum],[S.E90q],'b*');
%     p4=plot([S.Datenum],[S.E95q],'c*');
%     p5=plot([S.Datenum],[S.E98q],'r*');
%     p6=plot([S.Datenum],[S.Emax],'k^');
    for n=2:length([S.Datenum])
     p1=plot([S(n-1).Datenum S(n).Datenum],...
        [S(n).Emean S(n).Emean],'b-','linewidth',2);
     p2=plot([S(n-1).Datenum S(n).Datenum],...
        [S(n).Emedian S(n).Emedian],'g-','linewidth',2);
     p3=plot([S(n-1).Datenum S(n).Datenum],...
        [S(n).E90q S(n).E90q],'c-','linewidth',2);
     p4=plot([S(n-1).Datenum S(n).Datenum],...
        [S(n).E95q S(n).E95q],'y-','linewidth',2);
     p5=plot([S(n-1).Datenum S(n).Datenum],...
        [S(n).E98q S(n).E98q],'m-','linewidth',2);
     p6=plot([S(n-1).Datenum S(n).Datenum],...
        [S(n).Emax S(n).Emax],'r-','linewidth',2);
    end
    grid on;
    datetick(gca)
    ylabel('Wave Energy (m^2)');xlabel('Date');
    legend([p1 p2 p3 p4 p5 p6],'Mean','Median','Q90','Q95','Q98','Max',...
        'location','northwest');
    
    subplot(3,1,2)
    set(gca,'fontsize',14);hold on;
    plot([S(1).Datenum S(end).Datenum],[0 0],'k-');
    for n=2:length([S.Datenum])
     p1=plot([S(n-1).Datenum S(n).Datenum],...
        [S(n).SxyMean S(n).SxyMean],'b-','linewidth',2);
     p2=plot([S(n-1).Datenum S(n).Datenum],...
        [S(n).SxyMedian S(n).SxyMedian],'g-','linewidth',2);
     p3=plot([S(n-1).Datenum S(n).Datenum],...
        [S(n).SxyMax S(n).SxyMax],'r-','linewidth',2);
     plot([S(n-1).Datenum S(n).Datenum],...
        [S(n).SxyMin S(n).SxyMin],'r-','linewidth',2);
    end
    grid on;
    datetick(gca)
    ylabel('Sxy (m^2)');xlabel('Date');
    legend([p1 p2 p3],'Mean','Median','Min & Max',...
        'location','northwest');
      
    subplot(3,1,3)
    set(gca,'fontsize',14);hold on;
    plot([S(1).Datenum S(end).Datenum],[0 0],'k-');
    for n=2:length([S.Datenum])
     p1=plot([S(n-1).Datenum S(n).Datenum],...
        [S(n).SxyNet S(n).SxyNet],'b-','linewidth',2);
    end
    
    grid on;
    datetick(gca)
    ylabel('Cumulative Sxy (m^2)');xlabel('Date');
    %legend('Mean','Median','Max');

end 


fh=get(gcf);       

end


