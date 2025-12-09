load MopTableUTM.mat
load VosTransects2
% get unique mop numbers
vmops=round([VosT.BackMop]);
% sort unique mop numbers low to high
umops=sort(unique(vmops));

fid=fopen('CAmopSlopeParameters.txt','w');
S=46;

icount=0;
for n=umops
    B=median([VosT(vmops == n).Slope],'omitnan');
    % if no slope on related Vos transect, include a few neighbors
    if isnan(B) 
    fprintf('%5i %5s Null Slope, Searching Neigbors...\n',n,Mop.Name{n})
    B=median([VosT(vmops >= n-2 & vmops <= n+2).Slope],'omitnan');
    end
    % if still a Nan, skip
    if ~isnan(B)  
    A=B/2;
    fprintf(fid,'%5i %5s %5.3f %5.4f %3i\n',n,Mop.Name{n},B,A,S);
    icount=icount+1;
    else
    fprintf('%5i %5s No Nearby Slopes, skipping...\n',n,Mop.Name{n})   
    end
end
    
fclose(fid);
fprintf('Parameters for %i Mop Transects written to: %s\n',icount,'CAmopSlopeParameters.txt')
