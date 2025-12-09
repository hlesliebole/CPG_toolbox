
load /volumes/group/MOPS/M00582SG.mat

ndx=find(contains({SG.File},'umbo') == 1 );

fprintf('Index  Date      #grid_pts  minZ    maxZ\n')

for n=fliplr(ndx)
        fprintf('%4i %s %8i %6.1f %6.1f\n',n,datestr(SG(n).Datenum),...
        numel(SG(n).X),min(SG(n).Z),max(SG(n).Z))
end