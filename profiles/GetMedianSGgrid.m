function [X,Y,Zmed]=GetMedianSGgrid(SG)


%% returns the median grid matrix for the input SG struct array
%
%  X = 2d x utm indices
%  Y = 2d y utm indices
%  Zmed = median grid elevations

% medians are calculated are by progressivle claculating medians
%  based on year-month > month > quarter > season > global

%% Need to turn the saved SG struct array grid points into a 2d grid array

% make grid area encompassing all survey data
minx=min(vertcat(SG.X));
maxx=max(vertcat(SG.X));
miny=min(vertcat(SG.Y));
maxy=max(vertcat(SG.Y));
%[X,Y]=meshgrid(min(x):max(x),min(y):max(y));
[X,Y]=meshgrid(minx:maxx,miny:maxy);

Z=NaN(size(SG,2),size(X(:),1));

m=0;
for SurvNum=1:size(SG,2)
 m=m+1;
% Mop area 1m grid points with valid data
x=SG(SurvNum).X;
y=SG(SurvNum).Y;
% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);
Z1=X*NaN; % initialize the 2d elevation Z array as NaNs
% Z(m,idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices
% Zcount(idx)=Zcount(idx)+1;
% Zmean(idx)=Zmean(idx)+Z1(idx);
% Zmn=min(Zmn,Z1);
end

% survey datetimes
Zdt=datetime([SG.Datenum],'convertfrom','datenum');

% year-months
ym=year(Zdt)*100+month(Zdt);

% unique year-months
uym=unique(ym);

% loop through unique year-months and calculate a median
%  for each

m=0;
Z=NaN(numel(uym),size(X(:),1));
for i=uym
    m=m+1;
    idx=find(ym == i);
    j=0;
    Zym=NaN(numel(idx),size(X(:),1));
    for n=idx
      x=SG(SurvNum).X;
      y=SG(SurvNum).Y;
      % Make 2d x,y utm arrays encompassing the valid points
      ndx=sub2ind(size(X),y-miny+1,x-minx+1);
      j=j+1;
      Zym(j,idx)=SG(SurvNum).Z;
    end
      Z(m,:)=median(Zym,2);
end
    
    fprintf('%i\n',n)
end
    


% year-month means
n=0;
for y=year(Zdt(1)):year(Zdt(end))
   for m=1:12
       n=n+1;
       Zyear(n)=y;
       Zmonth(n)=m;
       Zym(n,:)=Zf(1,:)*NaN;
       idx=find(year(Zdt) == y & month(Zdt) == m);
       if numel(idx) == 1
           Zym(n,:)=Zf(idx,:);
       elseif numel(idx) > 1
           Zym(n,:)=mean(Zf(idx,:),'omitnan');
       end
   end
end

% monthly means across all years
n=0;
for m=1:12
       n=n+1;
       Zm(n,:)=Zym(1,:)*NaN;
       idx=find(Zmonth == m);
       if numel(idx) == 1
           Zm(n,:)=Zym(idx,:);
       elseif numel(idx) > 1
           Zm(n,:)=mean(Zym(idx,:),'omitnan');
       end
end


% quarterly means from monthly means
n=0;
for q=1:4
       n=n+1;
       Zq(n,:)=Zm(1,:)*NaN;      
       idx=(q-1)*3+1:q*3;
       Zq(n,:)=mean(Zm(idx,:),'omitnan');
end

% figure;
% hold on;
% for m=1:4
%     plot(X1Dt,Zq(m,:),'-','linewidth',2)
%     hold on;
% end
% legend

% seasonal means from quarterly means
Zs(1,:)=mean(Zq(1:2,:),'omitnan');
Zs(2,:)=mean(Zq(3:4,:),'omitnan');

% figure;
% hold on;
% for m=1:2
%     plot(X1Dt,Zs(m,:),'-','linewidth',2)
%     hold on;
% end
% legend
    
% global mean
Zg=mean(Zs,'omitnan');


end