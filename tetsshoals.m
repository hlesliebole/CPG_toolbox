z=zsave;
% find USACE Shoals survey indices
ShoalsIdx=find(strcmp({SG.Source},'USACE'));

% They are ordered oldest to newest so loop through
%  gridded shoals survey data and replace any 
%  base grid points with valid shoals gridded points. 

for sn=ShoalsIdx
  % 1m gridded survey x,y,z points with data
  xs=SG(sn).X;
  ys=SG(sn).Y;
  zs=SG(sn).Z;
  
  % gridded shoals data cab be weird around the surfzone
  %  and is mostly wanted for deeper depths so just use 
  %  data deeper than -5m
  
  idx=find( xs >= xmin & xs < xmax & ...
      ys >= ymin & ys < ymax & ...
      zs < -5);
  
  xs=xs(idx);
  ys=ys(idx);
  zs=zs(idx);
  
  % get 1d indices of xs,ys points in the base grid 
  idx=sub2ind(size(z),xs-xmin+1,ys-ymin+1);
  
  % replace these base grid data points
  z(idx)=zs; 

end

% make a figure of 1m spatial res base grid
figure('position',[ 27   231   537   543]);
imagesc(xmin:xmax,ymin:ymax,z');
demcmap(z);colorbar
set(gca,'ydir','normal','dataaspectratio',[1 1 1])
title('Base Grid +Shoals');xlabel('E UTM');ylabel('N UTM');

%% -----------------------------------------------------------------
%    STEP 4: Now overlay the Jan 24 2023 jumbo survey 
% ------------------------------------------------------------------

JumboIdx=find(contains({SG.File},'umbo') & ...
    [SG.Datenum] > datenum(2014,10,5) & ...
    [SG.Datenum] <= datenum(2023,1,24));

for sn=JumboIdx
% 1m gridded survey x,y,z points with data
  xs=SG(sn).X;
  ys=SG(sn).Y;
  zs=SG(sn).Z;
  
  idx=find( xs >= xmin & xs < xmax & ys >= ymin & ys < ymax );
  xs=xs(idx);
  ys=ys(idx);
  zs=zs(idx);
  
  % get 1d indices of xs,ys points in the base grid 
  idx=sub2ind(size(z),xs-xmin+1,ys-ymin+1);
  
  % replace these base grid data points
  z(idx)=zs; 
end

figure('position',[ 37   221   537   543]);
imagesc(xmin:xmax,ymin:ymax,z');
demcmap(z);colorbar
set(gca,'ydir','normal','dataaspectratio',[1 1 1])
title('Base Grid +Shoals + JUmbos');xlabel('E UTM');ylabel('N UTM');

