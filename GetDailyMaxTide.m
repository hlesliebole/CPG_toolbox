
load LaJollaMslTide20002022.mat 
TideDay=floor(datenum(Tide.Datenum)); % tide day
i=0;
for n=min(TideDay):max(TideDay)
 i=i+1;
 nh=find(TideDay == n);
 [DayMax(i),ih]=max(Tide.MSL(nh)+0.774); % daily max
end

