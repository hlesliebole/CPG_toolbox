clearvars
%load M00578SG.mat
SG=CombineSGdata('/Users/William/Desktop/MOPS/',575,598);
for n=1:size(SG,2)
    idx=find(SG(n).Z < -8 | SG(n).Z > 4);
    if ~isempty(idx)
     SG(n).X(idx)=[];
     SG(n).Y(idx)=[];
     SG(n).Z(idx)=[];
    end
end
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

for k=1:12
    

ikx=find(mons == k);
if ~isempty(ikx)
    
% unique year-months
uym=unique(ym(ikx));

% loop through unique year-months and calculate a median
%  for each

m=0;
Z=NaN(numel(uym),size(X(:),1));
for i=uym
    m=m+1;
    idx=find(ym == i & mons == k);
    %fprintf('%i  %i\n',i,numel(idx))
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
    

      
% now reduce the year-month medians for the month to a monthly median
    if m > 1
     Zmon(k,:)=median(Z,'omitnan');
    else
     Zmon(k,:)=Z(:);
    end
    
    Zcount(k)=m;

else
    Zcount(k)=0;
end

fprintf('Month/Num Surveys/Num Combined Year-Month Surveys: %i / %i / %i\n',...
    k,numel(ikx),Zcount(k))

end

% reduce monthly medians to quarterly medians
Zqrt=NaN(4,size(X(:),1));
for n=1:4
 Zqrt(n,:)=median(Zmon((n-1)*3+1:(n-1)*3+3,:),'omitnan');
end

% switch to calculating means at the seasonal and global levels

% reduce quarterly medians to seasonal means of the medians
Zssn=NaN(2,size(X(:),1));
Zssn(1,:)=mean(Zqrt([1 4],:),'omitnan');
Zssn(2,:)=mean(Zqrt([2 3],:),'omitnan');

% global mean of the seasonal means 
Zglo=X*NaN;Zglo(:)=mean(Zssn,'omitnan');

%figure;imagesc(Zglo);colormap(jet);colorbar
surf(X,Y,Zglo);shading flat;colorbar
demcmap(Zglo)
    
Xglo=X;Yglo=Y;
save TP578to595globalXYZ.mat Xglo Yglo Zglo

Zlast=X*NaN;
n=size(SG,2);
x=SG(n).X;
y=SG(n).Y;
% Make 2d x,y utm arrays encompassing the valid points
ndx=sub2ind(size(X),y-miny+1,x-minx+1);
      
Zlast(ndx)=SG(n).Z;

figure;surf(X,Y,Zlast-Zglo);shading flat;polarmap;colorbar
view(2)

hold on;

% load M00578SA.mat
% 
% figure
% hold on;
% 
%     for n=1:size(SA,2)
%       Zsa=X*NaN;
%       x=SA(n).X;
%       y=SA(n).Y;
%       % Make 2d x,y utm arrays encompassing the valid points
%       ndx=sub2ind(size(X),y-miny+1,x-minx+1);
%       j=j+1;
%       Zsa(ndx)=SA(n).Z;
%       dZ=Zsa-Zglo;
%       idx=find(~isnan(dZ(:)));
%       plot(Zglo(idx),dZ(idx),'.')
%     end
% 
% % dZ=Zsa-Zglo;
% % idx=find(~isnan(dZ(:)));
% % figure;plot(Zglo(idx),dZ(idx),'.')
