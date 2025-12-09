%
% matlab script to read binary .rbd files
% unformatted, direct access, 3 real*4 numbers per record
% (decimal lon, decimal lat, depth)
%

% make a menu ist of .rbd files

!ls -1 data/* > rbd_list

fid=fopen('rbd_list','r');

j=0;
all=[];
while 1
 line=fgetl(fid);
 if ~isstr(line), break, end;
 j=j+1;
eval(['opt' int2str(j) '=line;']);
filename=['opt' int2str(j)];
all = [ all  ',' filename];
end

status=fclose(fid);

ttl='Select file name ';
eval(['nf=menu(ttl' all ');']);

filename=['opt' int2str(nf)];
rbf=eval(filename);
%fid=fopen(rbf,'r');
%S=fread(fid,[3,inf],'float');
%status=fclose(fid);
S=load(rbf);
S=S';

% flip tahiti
%S(1,:)=-S(1,:);

if ns == 0;
lon=S(1,:)';
lat=S(2,:)';
dep=S(3,:)';
ns=max(size(lon));
plot_ttl= rbf ;
else
lon=[ lon ; S(1,:)'];
lat=[ lat ; S(2,:)'];
dep=[ dep ; S(3,:)'];
ns=max(size(lon));
plot_ttl=[ plot_ttl  ' &  '  rbf ];
end
clear S

