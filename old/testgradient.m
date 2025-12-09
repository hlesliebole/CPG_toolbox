clf
Z=[0 0 0 0 0 0 0 0 0 0 0;...
   0 0 0 0 0 0 0 0 0 0 0;...
   0 0 0 0 0 0 0 0 0 0 0;...
   0 0 0 1 1 1 1 1 0 0 0;...
   0 0 0 1 1 1 1 1 0 0 0;...
   0 0 0 1 1 1 1 1 0 0 0;...
   0 0 0 0 0 0 0 0 0 0 0;...
   0 0 0 0 0 0 0 0 0 0 0;...
   0 0 0 0 0 0 0 0 0 0 0];

dx1=diff(Z,1,2);dx1(:,end+1)=dx1(:,end); % x grid elev differences
dy1=diff(Z,1,1);dy1(end+1,:)=dy1(end,:); % y grid elevation differences
sign1=sign(cosd(Mop.Normal(Mopnum)-(270-atan2d(dy1,dx1))));
dx2=diff(fliplr(Z),1,2);dx2(:,end+1)=dx2(:,end);dx2=fliplr(dx2);
dy2=diff(flipud(Z),1,1);dy2(end+1,:)=dy2(end,:);dy2=flipud(dy2);
sign2=-sign(cosd(Mop.Normal(Mopnum)-(270-atan2d(dy2,dx2))));

signc=sign(sign1.*sqrt(dx1.^2+dy1.^2)+sign2.*sqrt(dx2.^2+dy2.^2));
stpx=(abs(dx1)+abs(dx2))/2;
stpy=(abs(dy1)+abs(dy2))/2;
stp=signc.*sqrt(stpx.^2+stpy.^2);

[fx,fy]=gradient(Z);
s=sqrt(fx.^2+fy.^2);
fmax=max(cat(3,abs(fx),abs(fy)),[],3);


% plot(y);hold on;
% plot(gradient(y),'--');
% plot(gradient(y),'x');
% d1=[diff(y) nan];
% plot(diff(y),'-.');
% d2=[nan fliplr(diff(fliplr(y)))];
% plot([nan fliplr(diff(fliplr(y)))],'k:');
% % plot(diff(gradient(y)),'o');
% %plot(2:1+length(diff(y,2)),diff(y,2));
% 
% dm=(abs(d1)+abs(d2))/2;  % mean steepness of the point
% plot(dm,'r-')