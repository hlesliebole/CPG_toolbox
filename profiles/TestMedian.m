clearvars
load M00582SG.mat
minx=min(vertcat(SG.X));
maxx=max(vertcat(SG.X));
miny=min(vertcat(SG.Y));
maxy=max(vertcat(SG.Y));
%[X,Y]=meshgrid(min(x):max(x),min(y):max(y));
[X,Y]=meshgrid(minx:maxx,miny:maxy);


Zdt=datetime([SG.Datenum],'convertfrom','datenum');

% initialize monthly median matrix
Zmon=NaN(12,size(X(:),1));

% year-months for this month
ym=year(Zdt)*100+month(Zdt);

% months
mons=month(Zdt);

% loop through unique months
umons=unique(mons);

for k=umons
    
fprintf('Month %i\n',k)
ikx=find(mons == k);
    
% unique year-months
uym=unique(ym(ikx));

% loop through unique year-months and calculate a median
%  for each

m=0;
Z=NaN(numel(uym),size(X(:),1));
for i=uym
    m=m+1;
    idx=find(ym == i & mons == k);
    fprintf('%i  %i\n',i,numel(idx))
    j=0;
    Zym=NaN(numel(idx),size(X(:),1));
    for n=idx
      x=SG(n).X;
      y=SG(n).Y;
      % Make 2d x,y utm arrays encompassing the valid points
      ndx=sub2ind(size(X),y-miny+1,x-minx+1);
      j=j+1;
      Zym(j,ndx)=SG(n).Z;
    end
      if j > 1
       Z(m,:)=median(Zym,'omitnan');
      else
       Z(m,:)=Zym;
      end
end
    
    Zmon(k,:)=median(Z,'omitnan');
    Zcount(k)=m;
    
end
    fprintf('%i\n',n)
    
    Zmed=X*NaN;Zmed(:)=median(Zmon,'omitnan');
    
