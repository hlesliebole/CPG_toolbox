prompt='VPN connection required; Enter Spot ID';
def={'?','hsv'};
idc=inputdlg(prompt,'Spot ID',1,def);
id=cell2mat(idc);
%id=input('Enter Spot ID:')
fid=fopen('spot.id','w');
fprintf(fid,'%g\n',id);
fclose(fid);
urlwrite(['http://www.surfline.com/xml/sh.cfm?id=' num2str(id)],'cms.dat');
!rbd_mfiles/get_cms_ll.sh;
fid=fopen('spot.dat','r');
spll=fscanf(fid,'%g %g\n');
fclose(fid)
plot(spll(1),spll(2),'m*','markersize',10);
plot(spll(1),spll(2),'ko','markersize',10);
