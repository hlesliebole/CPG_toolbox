load VosTransects.mat % load CoastSat transect info
for n=1:size(VosT,2)
    if isempty(VosT(n).BackMop) % fill any empty fields with NaNs
       VosT(n).BackMop=NaN;
       VosT(n).BackMopXoffset=NaN;
    end
end

idx=find([VosT.BackMop] > 0);

figure;plot([VosT(idx).BackMop],[VosT(idx).Slope],'.')
VosT(idx(1120))
VosT(idx(1121))