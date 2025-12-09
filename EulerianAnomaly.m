
% For an input global mean Mop profile elevation and Mop reach,
%  calculates the individual survey elevation anomalies at the xshore 
%  location(s) of the input elevation on each Mop's global mean profile.

%---------
% Settings
%---------
clear all
ElvMsl=2.5; % mean Mop profile elevation MSL for the anomaly estimates
Elv=ElvMsl+0.774; % convert to navd88

MopStart=520;  % mop range
MopEnd=596;
% MopStart=552;  % mop range
% MopEnd=552;
% running mean window (odd number of months) for mean seasonal MOP profiles
SeasonalAvgWindow=3;  

%----------------------

load MopTableUTM.mat

MopNumbers=MopStart:MopEnd;
NumMopPoints=length(MopNumbers);

GS.MopNumbers=MopNumbers;

% loop through Mop SM mat files

mm=0;
for m=MopStart:MopEnd  % Mop number loop
    mm=mm+1;
   
    % load mop SM mat file
    
    load([ 'M'  num2str( m , '%5.5i' )  'SM.mat' ],'SM');
    load([ 'M'  num2str( m , '%5.5i' )  'GM.mat' ],'GM');
    
    fprintf('Loaded: %s %i\n',[ 'M'  num2str( m , '%5.5i' )  'SM.mat' ],size(SM,2));
    
    % get global mean and seaosnal mean profile Elv locations for this mop
       clear gx gsx
       
    % global Elv location(s) 
        xl=[GM.X1D(1) GM.X1D(end)];  
        gx.x=intersections(xl,[Elv Elv],GM.X1D,GM.Z1Dmean);
        % round xshore locatiosn to nearest meter
        gx.x=round(gx.x);
     
     % seasonal mean Elv location(s)
          
          for mon=1:12
              % make mean seasonal profile from monthly mean
              % profiles using SeasonalAvgWindow setting
              dm=(SeasonalAvgWindow-1)/2;
              mons=mon-dm:mon+dm;
              mons(mons < 1)=mons(mons < 1)+12;
              mons(mons > 12)=mons(mons > 12)-12;
              Zsm=nanmean(reshape([GM.MM(mons).Z1Dmean],[length(GM.MM(1).Z1Dmean),SeasonalAvgWindow]),2);
              xl=[GM.MM(mon).X1D(1) GM.MM(mon).X1D(end)];
              %dist=intersections(xl,[Elv Elv],GM.MM(mon).X1D,Zsm);
              gsx.MM(mon).x=intersections(xl,[Elv Elv],GM.MM(mon).X1D,Zsm);
              gsx.MM(mon).x=round(gsx.MM(mon).x);
              if isempty(gsx.MM(mon).x)
                  gsx.MM(mon).x=NaN;
              end
          end
        
        
    % step through survey dates
    for n=1:size(SM,2)
        sdate=SM(n).Datenum;
        mon=month(datetime(datestr(sdate))); % survey month
        stype=SM(n).Source;
        
        % get anomaly relative to global mean profile
        idx=find(ismember(round(SM(n).X1D),gx.x));
        ganom=nanmean(SM(n).Z1Dmean(idx)-Elv);
        % get anomaly relative to global mean seasonal profile
        idx=find(ismember(round(SM(n).X1D),gsx.MM(mon).x));
        gsanom=nanmean(SM(n).Z1Dmean(idx)-Elv);
   
    
    %     % find index of exist Shoreline struct array entry for this date
    %  and survey source if it exists
        if m == MopStart
          idx=[]; 
        else
         idx=find([Anom.Datenum] == sdate & strcmpi({Anom.Source}, stype)==1);
        end
    
    % if new date, add new entry to the Shoreline struct array and
    %  populate shoreline vectors with with NaNs for all Mops
      
       if isempty(idx)
           if m == MopStart && n == 1
               nn=1;
           else
               nn=size(Anom,2)+1;
           end
        Anom(nn).Datenum=sdate;
        Anom(nn).Source=stype;
        Anom(nn).MopNumbers=MopStart:MopEnd;
        Anom(nn).Global=nan(1,NumMopPoints);
        Anom(nn).Seasonal=nan(1,NumMopPoints);
        
        Anom(nn).Global(mm)=ganom;
        Anom(nn).Seasonal(mm)=gsanom;
       else
        Anom(idx).Global(mm)=ganom;
        Anom(idx).Seasonal(mm)=gsanom;
       end
    end
end

figure('position',[415   244   469   404],'name','Mop Mean Contour Anomalies');
x=repmat([Anom.Datenum],NumMopPoints,1);x=x(:)'; % make survey date for every points
y=[Anom.MopNumbers];
z=[Anom.Seasonal];
amax=ceil(max(abs(z)));
amax=2;
zscaled = 1+64*(z+amax)/(2*amax);
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points
                                     
cm =flipud(polarmap(64));

scp=scatter3(x(idx), y(idx), z(idx), 15, cm(ceil(zscaled(idx)),:), 'filled');
scp.MarkerFaceAlpha = .9;
scp.MarkerEdgeAlpha = .9;
view(-10,80)
colormap(flipud(polarmap(64)))
cb=colorbar;cb.Label.String='Elevation Anomaly (m)';
set(gca,'clim',[-amax amax]);
set(gca,'xlim',[min(x) max(x)]);
set(gca,'zlim',[min(z) max(z)]);
set(gca,'color',[.7 .7 .7]);
datetick;
xlabel('Survey Date');ylabel('Mop Number');zlabel('Anomaly (m)');
title(['Elevation Anomalies at the ' num2str(ElvMsl) 'm MSL Seasonal Mean Contour']); 
view(2)
set(gcf,'inverthardcopy','off');

%         % Survey elevation and mean profile Elv location
%         xl=[SM(n).X1D(1) SM(n).X1D(end)];  
%         xMSL=intersections(xl,[-8 -8],SM(n).X1D,SM(n).Z1Dmean);
%         xMSL=max(xMSL);
%         if isempty(xMSL)
%             xMSL=NaN;
%             UTMxMSL=NaN;
%             UTMyMSL=NaN;
%             LatMSL=NaN;
%             LonMSL=NaN;
%         else
%             UTMxMSL=Mop.BackXutm(m)-xMSL*cos(MopTransectAngles(mm));
%             UTMyMSL=Mop.BackYutm(m)-xMSL*sin(MopTransectAngles(mm));
%             [LatMSL,LonMSL]=utm2deg(UTMxMSL,UTMyMSL,...
%             repmat(SM(n).UTMzone,[length(UTMxMSL) 1]));
%         end
%         
%         % MHW location
%         xMHW=intersections(xl,[1.344 1.344],SM(n).X1D,SM(n).Z1Dmean);
%         xMHW=max(xMHW);
%         if isempty(xMHW)
%             xMHW=NaN;
%             UTMxMHW=NaN;
%             UTMyMHW=NaN;
%             LatMHW=NaN;
%             LonMHW=NaN;
%         else
%             UTMxMHW=Mop.BackXutm(m)-xMHW*cos(MopTransectAngles(mm));
%             UTMyMHW=Mop.BackYutm(m)-xMHW*sin(MopTransectAngles(mm));
%             [LatMHW,LonMHW]=utm2deg(UTMxMHW,UTMyMHW,...
%             repmat(SM(n).UTMzone,[length(UTMxMHW) 1]));
%         end
%         
%     
%      % convert beach widths to UTM shoreline locations
%     
%     % find index of exist Shoreline struct array entry for this date
%     %  and survey source if it exists
%         if m == MopStart
%           idx=[]; 
%         else
%          idx=find([SS.Datenum] == sdate & strcmpi({SS.Source}, stype)==1);
%         end
%     
%     % if new date, add new entry to the Shoreline struct array and
%     %  populate shoreline vectors with with NaNs for all Mops
%        if isempty(idx)
%            if m == MopStart && n == 1
%                nn=1;
%            else
%                nn=size(SS,2)+1;
%            end
%            
%            % ad new entry
%            fprintf('Adding : %s %s\n',datestr(sdate),stype)
%            SS(nn).Datenum=sdate;
%            SS(nn).Source=stype;
%            SS(nn).UTMzone=SM(n).UTMzone;
%            SS(nn).MopNumbers=MopNumbers;
%            SS(nn).MSLlatitudes=nan(1,NumShorePoints);
%            SS(nn).MSLlongitudes=nan(1,NumShorePoints);
%            SS(nn).MSLeastings=nan(1,NumShorePoints);
%            SS(nn).MSLnorthings=nan(1,NumShorePoints);
%            SS(nn).MSLbeachWidths=nan(1,NumShorePoints);
%            SS(nn).MSLanomalies=nan(1,NumShorePoints);
%            SS(nn).MSLmonthanomalies=nan(1,NumShorePoints);
%            SS(nn).MHWlatitudes=nan(1,NumShorePoints);
%            SS(nn).MHWlongitudes=nan(1,NumShorePoints);
%            SS(nn).MHWeastings=nan(1,NumShorePoints);
%            SS(nn).MHWnorthings=nan(1,NumShorePoints);
%            SS(nn).MHWbeachWidths=nan(1,NumShorePoints);
%            SS(nn).MHWanomalies=nan(1,NumShorePoints);
%            SS(nn).MHWmonthanomalies=nan(1,NumShorePoints);
%            
%            SS(nn).MSLbeachWidths(mm)=xMSL;
%            SS(nn).MSLanomalies(mm)=xMSL-gxMSL;
%            mon=month(datetime(datestr(sdate)));
%            SS(nn).MSLmonthanomalies(mm)=xMSL-gmxMSL(mon);
%            SS(nn).MHWbeachWidths(mm)=xMHW;
%            SS(nn).MHWanomalies(mm)=xMHW-gxMHW;
%            SS(nn).MHWmonthanomalies(mm)=xMHW-gmxMHW(mon);
%            SS(nn).MSLeastings(mm)=UTMxMSL;
%            SS(nn).MSLnorthings(mm)=UTMyMSL;
%            SS(nn).MHWeastings(mm)=UTMxMHW;
%            SS(nn).MHWnorthings(mm)=UTMyMHW;
%            SS(nn).MSLlatitudes(mm)=LatMSL;
%            SS(nn).MSLlongitudes(mm)=LonMSL;
%            SS(nn).MHWlatitudes(mm)=LatMHW;
%            SS(nn).MHWlongitudes(mm)=LonMHW;
%            
%            
%        else
%            % if survey data already exists add shoreline info
%            SS(idx).MSLbeachWidths(mm)=xMSL;
%            SS(idx).MSLanomalies(mm)=xMSL-gxMSL;
%            mon=month(datetime(datestr(sdate)));
%            SS(idx).MSLmonthanomalies(mm)=xMSL-gmxMSL(mon);
%            SS(idx).MHWbeachWidths(mm)=xMHW;
%            SS(idx).MHWanomalies(mm)=xMHW-gxMHW;
%            SS(idx).MHWmonthanomalies(mm)=xMHW-gmxMHW(mon);
%            SS(idx).MSLeastings(mm)=UTMxMSL;
%            SS(idx).MSLnorthings(mm)=UTMyMSL;
%            SS(idx).MHWeastings(mm)=UTMxMHW;
%            SS(idx).MHWnorthings(mm)=UTMyMHW;
%            SS(idx).MSLlatitudes(mm)=LatMSL;
%            SS(idx).MSLlongitudes(mm)=LonMSL;
%            SS(idx).MHWlatitudes(mm)=LatMHW;
%            SS(idx).MHWlongitudes(mm)=LonMHW;  
%        end
%        
%     end %survey date loop
%     
% end % end Mop number loop
% 
% % sort the Shoreline struct array by survey date 
% T=struct2table(SS); % sort by date before saving
% sortedT = sortrows(T, 'Datenum');
% SS=table2struct(sortedT)';
% 
% % save in a mat file
% save TorreyPinesShoreline.mat SS GS
% 
% % MSL2d=reshape([SS.MSLanomalies],[length(SS(1).MSLbeachWidths),size(SS,2)]);
% % pcolor([SS.Datenum],MopNumbers,MSL2d);datetick;colormap(jet);colorbar;
% % xlabel('Date');ylabel('Mop Number');shading flat;
% % 
% % figure
% % for i=1:size(SS,2)
% %     %plot(SS(i).MHWeastings,SS(i).MHWnorthings,'-')
% %     plot(SS(i).MHWlongitudes,SS(i).MHWlatitudes,'-')
% %     hold on;
% % end
% % for m=MopStart:MopEnd
% %     %plot([Mop.BackXutm(m) Mop.OffXutm(m)],[Mop.BackYutm(m) Mop.OffYutm(m)],'k-');
% %     plot([Mop.BackLon(m) Mop.OffLon(m)],[Mop.BackLat(m) Mop.OffLat(m)],'k-'); 
% % end
% % %xlabel('UTM Easting');ylabel('UTM Northing');
% % xlabel('Longitude');ylabel('Latitude');
% %   
