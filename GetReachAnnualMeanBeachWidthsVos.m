
function [xmbw,mbw]=GetReachAnnualMeanBeachWidthsVos(MopRange)

load VosShoreline2.mat
MopShoreline=VosShoreline;
MopShoreline=renamevars(MopShoreline,"BackMopDist","MslDist")

x=[];y=[];zw=[];za=[];
byrs=[];
bwid=[];

mnum=MopShoreline.MopNum;
for m=MopRange(1):MopRange(2)%666:682;%620:635;%638:662;%666:682 % Cardiff SB reach
    idx=find(mnum == m); % find this mop's rows in table
    if ~isempty(idx)
    d=MopShoreline.DateNum(idx); %survey datnums
    
    % shoreline width
   %xSL=MopShoreline.MhwDist(idx); % MHW width option
   xSL=MopShoreline.MslDist(idx); % MSL width option (less valid data) 
    TF=isoutlier(xSL);xSL(TF)=NaN; % remove outliers

    % reduce to individual year-month means (ymmeans) 
    dtime=datetime(d,'convertfrom','datenum');
    y=year(dtime);mn=month(dtime);
    uy=unique(y); %unique years
    ymmean(1:numel(uy(1):uy(end)),1:12)=NaN; % initialize matrix as NaNs
    yrs=uy(1):uy(end); %loop through unique years
    for ny=uy'
        for m=unique(mn(y == ny))' % loop through unique months in year
            ymmean(ny-uy(1)+1,m)=mean(xSL(y == ny & mn == m),'omitnan');
        end
    end
   
    % use year-month means to calculate beach year annual means
    for y=yrs(2:end)
        iy=find(yrs == y); % find this year
        ipy=find(yrs == y-1); % find previous year
        if ~isempty(ipy)
            q1=mean(ymmean(ipy,9:12),'omitnan'); % (sep)Oct-Dec prev year
            q2=mean(ymmean(iy,1:3),'omitnan'); % Jan-Mar year
            q3=mean(ymmean(iy,4:6),'omitnan'); % Apr-Jun year
            q4=mean(ymmean(iy,7:10),'omitnan'); % Jul-Sep year
        end
            wmin=mean([q1 q2],'omitnan'); % winter mean
            smean=mean([q3 q4],'omitnan'); % summer mean
        % make an annual mean if have both winter and summer values    
        if ~isnan(mean([wmin smean],'omitnan'))
            ymean=mean([wmin smean],'omitnan');
            % add to single vectors of data for all mops in reach
            byrs=[byrs y]; % beach year
            bwid=[bwid ymean]; % width
        end
    end
end
  
% Finally, spatially average means over the entire mop range
    fprintf('BeachYr   Width  #Mops with Data\n')
    xmbw=[];mbw=[];
    uy=unique(byrs);
    for y=uy % unique beach year loop
        iy=find(byrs == y);
        mbwid=mean(bwid(iy),'omitnan'); % spatial average
        fprintf('%i     %6.2f        %i\n',y,mbwid,numel(iy))
        % only return a value for a year if > 75% of mop reach has valid widths
        if numel(iy)/(MopRange(2)-MopRange(1)) > .75
        xmbw=[xmbw y];
        mbw=[mbw mbwid];
        end
    end
end
end
    
