function [uy,ymmeans]=GetSAtimeHistory(SA)



dtime=datetime([SA.Datenum],'convertfrom','datenum');
y=year(dtime);
mn=month(dtime);
dy=day(dtime);

fprintf('%s\n',...
'  YEAR   Jan   Feb   Mar   Apr   May   Jun   Jul   Aug   Sep   Oct   Nov   Dec');
% Reduce to counts of individual year-month means
uy=unique(y); %unique years
ymmeans(1:numel(uy(1):uy(end)),1:12)=0;
for ny=uy
    for m=unique(mn(y == ny))
        ymmeans(ny-uy(1)+1,m)=numel(find(y == ny & mn == m));
    end
    fprintf('%6i',[ny ymmeans(ny-uy(1)+1,:)]')
    fprintf('\n')
end
        
end