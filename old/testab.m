
m=513;
load MopTableUTM.mat
load(['M' num2str(m,'%5.5i') 'SA.mat']);

a=(Mop.BackYutm(m)-Mop.OffYutm(m))/...
    (Mop.BackXutm(m)-Mop.OffXutm(m));

b=Mop.BackYutm(m)-a*Mop.BackXutm(m);

d=DistanceFromLine(SA(44).X,SA(44).Y,a,b);
figure;plot(SA(44).Z,d,'k.')

k = ((y2-y1) * (x3-x1) - (x2-x1) * (y3-y1)) / ((y2-y1)^2 + (x2-x1)^2)
x4 = x3 - k * (y2-y1)
y4 = y3 + k * (x2-x

function d=DistanceFromLine(x,y,a,b)

% perpendicular distance of the point x,y from the line defined
%   by y = ax + b

d=abs(a.*x - y + b)./sqrt(a.^2+1);

end