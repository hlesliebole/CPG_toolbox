%
% matlab script to prompt and load a M*SA.mat
%  file and its accompanying M*QC.mat file if it exists

% make a figure with Mop number input box in upper left corner



% default start is 654 fletcher cove
%CurrentMopNumber=654;
% MopNumber=654;
% load('M00654SA.mat');
axes(ax1)

if exist('nav','var');delete(nav);end


Mhead=uicontrol(MainFig,'style','text','position',[10 660 70 30],...
    'string','MOP #','foregroundcolor','b','backgroundcolor','w',...
    'fontsize',22);

Medit=uicontrol(MainFig,'style','edit','position',[85 660 100 30],...
    'string',CurrentMopNumber,'Value',1,'fontsize',22,'backgroundcolor',[.9 .9 .9],...
    'callback',...
    'CurrentMopNumber=str2num(Medit.String);delete(Mhead);delete(Medit);GetMopMobileLidarSA');


% 
% !ls -1 data/* > rbd_list
% 
% fid=fopen('rbd_list','r');
% 
% j=0;
% all=[];
% while 1
%  line=fgetl(fid);
%  if ~isstr(line), break, end;
%  j=j+1;
% eval(['opt' int2str(j) '=line;']);
% filename=['opt' int2str(j)];
% all = [ all  ',' filename];
% end
% 
% status=fclose(fid);
% 
% ttl='Select file name ';
% eval(['nf=menu(ttl' all ');']);
% 
% filename=['opt' int2str(nf)];
% rbf=eval(filename);
% %fid=fopen(rbf,'r');
% %S=fread(fid,[3,inf],'float');
% %status=fclose(fid);
% S=load(rbf);
% S=S';
% 
% % flip tahiti
% %S(1,:)=-S(1,:);
% 
% if ns == 0;
% lon=S(1,:)';
% lat=S(2,:)';
% dep=S(3,:)';
% ns=max(size(lon));
% plot_ttl= rbf ;
% else
% lon=[ lon ; S(1,:)'];
% lat=[ lat ; S(2,:)'];
% dep=[ dep ; S(3,:)'];
% ns=max(size(lon));
% plot_ttl=[ plot_ttl  ' &  '  rbf ];
% end
% clear S
% 
