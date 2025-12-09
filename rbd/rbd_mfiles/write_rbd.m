%
% matlab script to write binary edited.rbd file
% unformatted, direct access, 3 real*4 numbers per record
% (decimal lon, decimal lat, depth)
%

% flip sign of longitude for storing   

save temp.mat lon lat    

S(1,:)= lon';
clear lon
S(2,:)= lat';
clear lat
S(3,:)= dep';

fid=fopen('data/EDITED.xyz','w');
%w=fwrite(fid,S,'float');
fprintf(fid,'%10.4f %10.4f %8.1f\n',S);
fclose(fid);

clear S;

load temp
%alon=lon;
%alat=lat;
%adep=dep;

%clear alon alat adep

'Edited data written to data/EDITED.xyz '
