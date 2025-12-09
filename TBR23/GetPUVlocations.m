% get lat/lons of -5 and -7 msl locations on mops
% 580 and 586

load M00580SA.mat
idx=find([SA.Datenum] == datenum(2023,4,28));

fprintf('\n     Mop 580\n\n')
% -7m msl
ndx=find(SA(idx).Z < (-7+0.774));
[zmin,imax]=max(SA(idx).Z(ndx));
SA(idx).Z(ndx(imax));
[lat7,lon7]=utm2deg(SA(idx).X(ndx(imax)),SA(idx).Y(ndx(imax)),'11 S');
fprintf('7m MSL: %8.5f N %12.5f W \n',lat7,lon7)

% -5m msl
ndx=find(SA(idx).Z < (-5+0.774));
[zmin,imax]=max(SA(idx).Z(ndx));
SA(idx).Z(ndx(imax));
[lat5,lon5]=utm2deg(SA(idx).X(ndx(imax)),SA(idx).Y(ndx(imax)),'11 S');
fprintf('5m MSL: %8.5f N %12.5f W \n',lat5,lon5)

load M00586SA.mat
idx=find([SA.Datenum] == datenum(2023,4,28));

fprintf('\n     Mop 586\n\n')
% -7m msl
ndx=find(SA(idx).Z < (-7+0.774));
[zmin,imax]=max(SA(idx).Z(ndx));
SA(idx).Z(ndx(imax));
SA(idx).Y(ndx(imax));
[lat7,lon7]=utm2deg(SA(idx).X(ndx(imax)),SA(idx).Y(ndx(imax)),'11 S');
fprintf('7m MSL: %8.5f N %12.5f W \n',lat7,lon7)

% -5m msl
ndx=find(SA(idx).Z < (-5+0.774));
[zmin,imax]=max(SA(idx).Z(ndx));
SA(idx).Z(ndx(imax));
SA(idx).Y(ndx(imax));
[lat5,lon5]=utm2deg(SA(idx).X(ndx(imax)),SA(idx).Y(ndx(imax)),'11 S');
fprintf('5m MSL: %8.5f N %12.5f W \n',lat5,lon5)