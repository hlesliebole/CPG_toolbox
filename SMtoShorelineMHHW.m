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
% SS.MHHWlatitudes
% SS.MHHWlongitudes
% SS.MHHWeastings
% SS.MHHWnorthings
% SS.MHHWbeachWidths

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
        
        
        % MHHW location
        gxMHHW=intersections(xl,[1.566 1.566],GM.X1D,GM.Z1Dmean);
        gxMHHW=max(gxMHHW);
        if isempty(gxMHHW)
            gxMHHW=NaN;
            gxVolume=NaN;
            gxCentroid=NaN;
            gxSlope=NaN;
            gUTMxMHHW=NaN;
            gUTMyMHHW=NaN;
            gLatMHHW=NaN;
            gLonMHHW=NaN;
        else
            if isempty(gxMSL)
                gxVolume=NaN;
                gxCentroid=NaN;
                gxSlope=NaN;
            end
            [gxVolume,gxCentroid,gxSlope]=GetShoreface(gxMSL,gxMHHW,GM.X1D,GM.Z1Dmean);
            gUTMxMHHW=Mop.BackXutm(m)-gxMHHW*cos(MopTransectAngles(mm));
            gUTMyMHHW=Mop.BackYutm(m)-gxMHHW*sin(MopTransectAngles(mm));
            [gLatMHHW,gLonMHHW]=utm2deg(gUTMxMHHW,gUTMyMHHW,...
            repmat(SM(1).UTMzone,[length(gUTMxMHHW) 1]));
        end
        % also get global month (seasonal) mean MHHW locations
          for mon=1:12
              xl=[GM.MM(mon).X1D(1) GM.MM(mon).X1D(end)];
              dist=intersections(xl,[1.344 1.344],...
                  GM.MM(mon).X1D,GM.MM(mon).Z1Dmean);
              if isempty(dist)
                  gmxMHHW(mon)=NaN;
                  gmxVolume(mon)=NaN;
                  gmxCentroid(mon)=NaN;
                  gmxSlope(mon)=NaN;
              else
                  gmxMHHW(mon)=max(dist);
                  if isempty(gmxMSL(mon))
                    gmxVolume(mon)=NaN;
                    gmxCentroid(mon)=NaN;
                    gmxSlope(mon)=NaN;
                  else
                  [gmxVolume(mon),gmxCentroid(mon),gmxSlope(mon),SlopeLF,rmseLF]...
                      =GetShoreface(gmxMSL(mon),gmxMHHW(mon),GM.MM(mon).X1D,GM.MM(mon).Z1Dmean);
                  end
              end
        
          end
        
           GS.MSLbeachWidths(mm)=gxMSL;
           GS.MHHWbeachWidths(mm)=gxMHHW;
           GS.MSLeastings(mm)=gUTMxMSL;
           GS.MSLnorthings(mm)=gUTMyMSL;
           GS.MHHWeastings(mm)=gUTMxMHHW;
           GS.MHHWnorthings(mm)=gUTMyMHHW;
           GS.MSLlatitudes(mm)=gLatMSL;
           GS.MSLlongitudes(mm)=gLonMSL;
           GS.MHHWlatitudes(mm)=gLatMHHW;
           GS.MHHWlongitudes(mm)=gLonMHHW;
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
        
        % MHHW location
        if numel(SM(n).Z1Dmean) > 3
        xMHHW=intersections(xl,[1.566 1.566],SM(n).X1D,SM(n).Z1Dmean);
        xMHHW=max(xMHHW);
        else
            xMHHW=[];
            
        end
        if isempty(xMHHW)
            xMHHW=NaN;
            Centroid=NaN;
            Volume=NaN;
            Slope=NaN;
            UTMxMHHW=NaN;
            UTMyMHHW=NaN;
            LatMHHW=NaN;
            LonMHHW=NaN;
        else
            if isempty(xMSL)
             Centroid=NaN;
             Volume=NaN;
             Slope=NaN;
            else
             [Volume,Centroid,Slope,SLopeLF,rmseLF]=GetShoreface(xMSL,xMHHW,SM(n).X1D,SM(n).Z1Dmean);
            end
            UTMxMHHW=Mop.BackXutm(m)-xMHHW*cos(MopTransectAngles(mm));
            UTMyMHHW=Mop.BackYutm(m)-xMHHW*sin(MopTransectAngles(mm));
            [LatMHHW,LonMHHW]=utm2deg(UTMxMHHW,UTMyMHHW,...
            repmat(SM(n).UTMzone,[length(UTMxMHHW) 1]));
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
           SS(nn).MHHWlatitudes=nan(1,NumShorePoints);
           SS(nn).MHHWlongitudes=nan(1,NumShorePoints);
           SS(nn).MHHWeastings=nan(1,NumShorePoints);
           SS(nn).MHHWnorthings=nan(1,NumShorePoints);
           SS(nn).MHHWbeachWidths=nan(1,NumShorePoints);
           SS(nn).MHHWanomalies=nan(1,NumShorePoints);
           SS(nn).MHHWmonthanomalies=nan(1,NumShorePoints);
           SS(nn).Volume=nan(1,NumShorePoints);
           SS(nn).Centroid=nan(1,NumShorePoints);
           SS(nn).Slope=nan(1,NumShorePoints);
           
           SS(nn).MSLbeachWidths(mm)=xMSL;
           SS(nn).MSLanomalies(mm)=xMSL-gxMSL;
           mon=month(datetime(datestr(sdate)));
           SS(nn).MSLmonthanomalies(mm)=xMSL-gmxMSL(mon);
           SS(nn).MHHWbeachWidths(mm)=xMHHW;
           SS(nn).MHHWanomalies(mm)=xMHHW-gxMHHW;
           SS(nn).MHHWmonthanomalies(mm)=xMHHW-gmxMHHW(mon);
           SS(nn).MSLeastings(mm)=UTMxMSL;
           SS(nn).MSLnorthings(mm)=UTMyMSL;
           SS(nn).MHHWeastings(mm)=UTMxMHHW;
           SS(nn).MHHWnorthings(mm)=UTMyMHHW;
           SS(nn).MSLlatitudes(mm)=LatMSL;
           SS(nn).MSLlongitudes(mm)=LonMSL;
           SS(nn).MHHWlatitudes(mm)=LatMHHW;
           SS(nn).MHHWlongitudes(mm)=LonMHHW;
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
           SS(idx).MHHWbeachWidths(mm)=xMHHW;
           SS(idx).MHHWanomalies(mm)=xMHHW-gxMHHW;
           SS(idx).MHHWmonthanomalies(mm)=xMHHW-gmxMHHW(mon);
           SS(idx).MSLeastings(mm)=UTMxMSL;
           SS(idx).MSLnorthings(mm)=UTMyMSL;
           SS(idx).MHHWeastings(mm)=UTMxMHHW;
           SS(idx).MHHWnorthings(mm)=UTMyMHHW;
           SS(idx).MSLlatitudes(mm)=LatMSL;
           SS(idx).MSLlongitudes(mm)=LonMSL;
           SS(idx).MHHWlatitudes(mm)=LatMHHW;
           SS(idx).MHHWlongitudes(mm)=LonMHHW;
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
%     %plot(SS(i).MHHWeastings,SS(i).MHHWnorthings,'-')
%     plot(SS(i).MHHWlongitudes,SS(i).MHHWlatitudes,'-')
%     hold on;
% end
% for m=MopStart:MopEnd
%     %plot([Mop.BackXutm(m) Mop.OffXutm(m)],[Mop.BackYutm(m) Mop.OffYutm(m)],'k-');
%     plot([Mop.BackLon(m) Mop.OffLon(m)],[Mop.BackLat(m) Mop.OffLat(m)],'k-'); 
% end
% %xlabel('UTM Easting');ylabel('UTM Northing');
% xlabel('Longitude');ylabel('Latitude');
%   
