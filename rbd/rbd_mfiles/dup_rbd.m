%
% matlab script to 
% remove duplicate points by comparing adjacent values of lat and lon
% repeat until the array size stays the same after duplicate removal

' removing duplicate points '

ds=max(size(lon))
% make two new arrays that are shifted by on value from original
x=lon; x(1)=[]; x(max(size(x))+1)=0;
y=lat; y(1)=[]; y(max(size(y))+1)=0;
% calculate diff between shifted and nonshifted arrays
% z=0= duplicate points
z=abs(x-lon)+abs(y-lat);
% find retains indices of all nonzero (nonduplicate) points
I=find(z);
% wax duplicates
lon=lon(I);
lat=lat(I);
dep=dep(I);
% run again if duplicate points were found (in case of triplicates
%  etc.)
de=max(size(lon))
% it should also catch triplicates etc.
%if de < ds, dup_rbd, end

clear I x y z ;
