function SS=SMtoShoreline(MopStart,MopEnd)
%MopStart=580;MopEnd=584;
% SMtoShoreline.m

% Makes a struct array of Mop transect shoreline location info
% for a specified range of Mops

% SS.Datenum
% SS.Source
% SS.UTMzone
% SS.MSLlatitudes
% SS.MSLlongitudes
% SS.MSLeastings
% SS.MSLnorthings
% SS.MSLbeachWidths
% SS.MHWlatitudes
% SS.MHWlongitudes
% SS.MHWeastings
% SS.MHWnorthings
% SS.MHWbeachWidths

% close all
% clear all

load MopTableUTM.mat

% MopStart=520;
% MopEnd=596;

MopNumbers=MopStart:MopEnd;
MopTransectAngles=atan((Mop.BackYutm(MopNumbers)-Mop.OffYutm(MopNumbers))...
    ./(Mop.BackXutm(MopNumbers)-Mop.OffXutm(MopNumbers)));
% number of shoreline points in the Shoreline struct array
NumShorePoints=length(MopNumbers);

GS.MopNumbers=MopNumbers;

% loop through Mop SM mat files

mm=0;
for m=MopStart:MopEnd  % Mop number loop
    mm=mm+1;
   
    % load mop SM mat file
    
    load([ 'M'  num2str( m , '%5.5i' )  'SM.mat' ],'SM');
    load([ 'M'  num2str( m , '%5.5i' )  'GM.mat' ],'GM');
    
    fprintf('Loaded: %s %i\n',[ 'M'  num2str( m , '%5.5i' )  'SM.mat' ],size(SM,2));
    
    % get global mean profiles shoreline locations for this mop
    
    % global MSL location 
        xl=[GM.X1D(1) GM.X1D(end)];  
        gxMSL=intersections(xl,[0.774 0.774],GM.X1D,GM.Z1Dmean);
        gxMSL=max(gxMSL);
        if isempty(gxMSL)
            gxMSL=NaN;
            gUTMxMSL=NaN;
            gUTMyMSL=NaN;
            gLatMSL=NaN;
            gLonMSL=NaN;
        else
            gUTMxMSL=Mop.BackXutm(m)-gxMSL*cos(MopTransectAngles(mm));
            gUTMyMSL=Mop.BackYutm(m)-gxMSL*sin(MopTransectAngles(mm));
            [gLatMSL,gLonMSL]=utm2deg(gUTMxMSL,gUTMyMSL,...
            repmat(SM(1).UTMzone,[length(gUTMxMSL) 1]));
        end
        % also get global month (seasonal) mean MSL locations
          for mon=1:12
              xl=[GM.MM(mon).X1D(1) GM.MM(mon).X1D(end)];
              dist=intersections(xl,[0.774 0.774],...
                  GM.MM(mon).X1D,GM.MM(mon).Z1Dmean);
              if isempty(dist)
                  gmxMSL(mon)=NaN;
              else
                  gmxMSL(mon)=max(dist);
              end
          end
        
        
        % MHW location
        gxMHW=intersections(xl,[1.344 1.344],GM.X1D,GM.Z1Dmean);
        gxMHW=max(gxMHW);
        if isempty(gxMHW)
            gxMHW=NaN;
            gxVolume=NaN;
            gxCentroid=NaN;
            gxSlope=NaN;
            gUTMxMHW=NaN;
            gUTMyMHW=NaN;
            gLatMHW=NaN;
            gLonMHW=NaN;
        else
            if isempty(gxMSL)
                gxVolume=NaN;
                gxCentroid=NaN;
                gxSlope=NaN;
            end
            [gxVolume,gxCentroid,gxSlope]=GetShoreface(gxMSL,gxMHW,GM.X1D,GM.Z1Dmean);
            gUTMxMHW=Mop.BackXutm(m)-gxMHW*cos(MopTransectAngles(mm));
            gUTMyMHW=Mop.BackYutm(m)-gxMHW*sin(MopTransectAngles(mm));
            [gLatMHW,gLonMHW]=utm2deg(gUTMxMHW,gUTMyMHW,...
            repmat(SM(1).UTMzone,[length(gUTMxMHW) 1]));
        end
        % also get global month (seasonal) mean MHW locations
          for mon=1:12
              xl=[GM.MM(mon).X1D(1) GM.MM(mon).X1D(end)];
              dist=intersections(xl,[1.344 1.344],...
                  GM.MM(mon).X1D,GM.MM(mon).Z1Dmean);
              if isempty(dist)
                  gmxMHW(mon)=NaN;
                  gmxVolume(mon)=NaN;
                  gmxCentroid(mon)=NaN;
                  gmxSlope(mon)=NaN;
              else
                  gmxMHW(mon)=max(dist);
                  if isempty(gmxMSL(mon))
                    gmxVolume(mon)=NaN;
                    gmxCentroid(mon)=NaN;
                    gmxSlope(mon)=NaN;
                  else
                  [gmxVolume(mon),gmxCentroid(mon),gmxSlope(mon),SlopeLF,rmseLF]...
                      =GetShoreface(gmxMSL(mon),gmxMHW(mon),GM.MM(mon).X1D,GM.MM(mon).Z1Dmean);
                  end
              end
        
          end
        
           GS.MSLbeachWidths(mm)=gxMSL;
           GS.MHWbeachWidths(mm)=gxMHW;
           GS.MSLeastings(mm)=gUTMxMSL;
           GS.MSLnorthings(mm)=gUTMyMSL;
           GS.MHWeastings(mm)=gUTMxMHW;
           GS.MHWnorthings(mm)=gUTMyMHW;
           GS.MSLlatitudes(mm)=gLatMSL;
           GS.MSLlongitudes(mm)=gLonMSL;
           GS.MHWlatitudes(mm)=gLatMHW;
           GS.MHWlongitudes(mm)=gLonMHW;
           GS.Volume(mm)=gxVolume;
           GS.Centroid(mm)=gxCentroid;
           GS.Slope(mm)=gxSlope;
    
    % step through survey dates
    for n=1:size(SM,2)
        
        if isempty(SM(n).UTMzone)
        SM(n).UTMzone='11 S';
        end
        sdate=SM(n).Datenum;
        stype=SM(n).Source;
        
        % MSL location
        if numel(SM(n).Z1Dmean) > 3
        xl=[SM(n).X1D(1) SM(n).X1D(end)];  
        xMSL=intersections(xl,[0.774 0.774],SM(n).X1D,SM(n).Z1Dmean);
        xMSL=max(xMSL);
        else
            xMSL=[];
            SM(n).UTMzone='11 S';
        end
        if isempty(xMSL)
            xMSL=NaN;
            UTMxMSL=NaN;
            UTMyMSL=NaN;
            LatMSL=NaN;
            LonMSL=NaN;
        else
            UTMxMSL=Mop.BackXutm(m)-xMSL*cos(MopTransectAngles(mm));
            UTMyMSL=Mop.BackYutm(m)-xMSL*sin(MopTransectAngles(mm));
            [LatMSL,LonMSL]=utm2deg(UTMxMSL,UTMyMSL,...
            repmat(SM(n).UTMzone,[length(UTMxMSL) 1]));
        end
        
        % MHW location
        if numel(SM(n).Z1Dmean) > 3
        xMHW=intersections(xl,[1.344 1.344],SM(n).X1D,SM(n).Z1Dmean);
        xMHW=max(xMHW);
        else
            xMHW=[];
            
        end
        if isempty(xMHW)
            xMHW=NaN;
            Centroid=NaN;
            Volume=NaN;
            Slope=NaN;
            UTMxMHW=NaN;
            UTMyMHW=NaN;
            LatMHW=NaN;
            LonMHW=NaN;
        else
            if isempty(xMSL)
             Centroid=NaN;
             Volume=NaN;
             Slope=NaN;
            else
             [Volume,Centroid,Slope,SLopeLF,rmseLF]=GetShoreface(xMSL,xMHW,SM(n).X1D,SM(n).Z1Dmean);
            end
            UTMxMHW=Mop.BackXutm(m)-xMHW*cos(MopTransectAngles(mm));
            UTMyMHW=Mop.BackYutm(m)-xMHW*sin(MopTransectAngles(mm));
            [LatMHW,LonMHW]=utm2deg(UTMxMHW,UTMyMHW,...
            repmat(SM(n).UTMzone,[length(UTMxMHW) 1]));
        end
        
    
     % convert beach widths to UTM shoreline locations
    
    % find index of exist Shoreline struct array entry for this date
    %  and survey source if it exists
        if m == MopStart
          idx=[]; 
        else
         idx=find([SS.Datenum] == sdate & strcmpi({SS.Source}, stype)==1);
        end
    
    % if new date, add new entry to the Shoreline struct array and
    %  populate shoreline vectors with with NaNs for all Mops
       if isempty(idx)
           if m == MopStart && n == 1
               nn=1;
           else
               nn=size(SS,2)+1;
           end
           
           % ad new entry
           fprintf('Adding : %i %s %s\n',n,datestr(sdate),stype)
           SS(nn).Datenum=sdate;
           SS(nn).Source=stype;
           SS(nn).UTMzone=SM(n).UTMzone;
           SS(nn).MopNumbers=MopNumbers;
           SS(nn).MSLlatitudes=nan(1,NumShorePoints);
           SS(nn).MSLlongitudes=nan(1,NumShorePoints);
           SS(nn).MSLeastings=nan(1,NumShorePoints);
           SS(nn).MSLnorthings=nan(1,NumShorePoints);
           SS(nn).MSLbeachWidths=nan(1,NumShorePoints);
           SS(nn).MSLanomalies=nan(1,NumShorePoints);
           SS(nn).MSLmonthanomalies=nan(1,NumShorePoints);
           SS(nn).MHWlatitudes=nan(1,NumShorePoints);
           SS(nn).MHWlongitudes=nan(1,NumShorePoints);
           SS(nn).MHWeastings=nan(1,NumShorePoints);
           SS(nn).MHWnorthings=nan(1,NumShorePoints);
           SS(nn).MHWbeachWidths=nan(1,NumShorePoints);
           SS(nn).MHWanomalies=nan(1,NumShorePoints);
           SS(nn).MHWmonthanomalies=nan(1,NumShorePoints);
           SS(nn).Volume=nan(1,NumShorePoints);
           SS(nn).Centroid=nan(1,NumShorePoints);
           SS(nn).Slope=nan(1,NumShorePoints);
           
           SS(nn).MSLbeachWidths(mm)=xMSL;
           SS(nn).MSLanomalies(mm)=xMSL-gxMSL;
           mon=month(datetime(datestr(sdate)));
           SS(nn).MSLmonthanomalies(mm)=xMSL-gmxMSL(mon);
           SS(nn).MHWbeachWidths(mm)=xMHW;
           SS(nn).MHWanomalies(mm)=xMHW-gxMHW;
           SS(nn).MHWmonthanomalies(mm)=xMHW-gmxMHW(mon);
           SS(nn).MSLeastings(mm)=UTMxMSL;
           SS(nn).MSLnorthings(mm)=UTMyMSL;
           SS(nn).MHWeastings(mm)=UTMxMHW;
           SS(nn).MHWnorthings(mm)=UTMyMHW;
           SS(nn).MSLlatitudes(mm)=LatMSL;
           SS(nn).MSLlongitudes(mm)=LonMSL;
           SS(nn).MHWlatitudes(mm)=LatMHW;
           SS(nn).MHWlongitudes(mm)=LonMHW;
           SS(nn).Volume(mm)=Volume;
           SS(nn).Centroid(mm)=Centroid;
           SS(nn).Slope(mm)=Slope;
           
           
       else
           idx=idx(end);
           % if survey data already exists add shoreline info
           SS(idx).MSLbeachWidths(mm)=xMSL;
           SS(idx).MSLanomalies(mm)=xMSL-gxMSL;
           mon=month(datetime(datestr(sdate)));
           SS(idx).MSLmonthanomalies(mm)=xMSL-gmxMSL(mon);
           SS(idx).MHWbeachWidths(mm)=xMHW;
           SS(idx).MHWanomalies(mm)=xMHW-gxMHW;
           SS(idx).MHWmonthanomalies(mm)=xMHW-gmxMHW(mon);
           SS(idx).MSLeastings(mm)=UTMxMSL;
           SS(idx).MSLnorthings(mm)=UTMyMSL;
           SS(idx).MHWeastings(mm)=UTMxMHW;
           SS(idx).MHWnorthings(mm)=UTMyMHW;
           SS(idx).MSLlatitudes(mm)=LatMSL;
           SS(idx).MSLlongitudes(mm)=LonMSL;
           SS(idx).MHWlatitudes(mm)=LatMHW;
           SS(idx).MHWlongitudes(mm)=LonMHW;
           SS(nn).Volume(mm)=Volume;
           SS(nn).Centroid(mm)=Centroid;
           SS(nn).Slope(mm)=Slope;
       end
       
    end %survey date loop
    
end % end Mop number loop

% sort the Shoreline struct array by survey date 
T=struct2table(SS); % sort by date before saving
sortedT = sortrows(T, 'Datenum');
SS=table2struct(sortedT)';

% save in a mat file
%save TorreyPinesShoreline.mat SS GS

end
% MSL2d=reshape([SS.MSLanomalies],[length(SS(1).MSLbeachWidths),size(SS,2)]);
% pcolor([SS.Datenum],MopNumbers,MSL2d);datetick;colormap(jet);colorbar;
% xlabel('Date');ylabel('Mop Number');shading flat;
% 
% figure
% for i=1:size(SS,2)
%     %plot(SS(i).MHWeastings,SS(i).MHWnorthings,'-')
%     plot(SS(i).MHWlongitudes,SS(i).MHWlatitudes,'-')
%     hold on;
% end
% for m=MopStart:MopEnd
%     %plot([Mop.BackXutm(m) Mop.OffXutm(m)],[Mop.BackYutm(m) Mop.OffYutm(m)],'k-');
%     plot([Mop.BackLon(m) Mop.OffLon(m)],[Mop.BackLat(m) Mop.OffLat(m)],'k-'); 
% end
% %xlabel('UTM Easting');ylabel('UTM Northing');
% xlabel('Longitude');ylabel('Latitude');
%   
