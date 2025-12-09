function [Volume,Centroid,Slope,SlopeLF,rmseLF]=GetShoreface(x1,x2,X1D,Z1Dmean)
  xq=round(X1D(1),2):0.01:round(X1D(end),2);
  zq=interp1(X1D,Z1Dmean,xq);
  idx=find(xq >= min([x1 x2]) & xq <= max([x1 x2]));
  Volume=trapz(xq(idx),zq(idx));
  Centroid=trapz(zq(idx).*xq(idx))/Volume;
  Slope=(1.566-0.774)/abs(x1-x2);
  [a,b]= linfit(xq(idx)',zq(idx)');
  SlopeLF=-a;
  rmseLF=sqrt(var(zq(idx)-(a.*xq(idx)+b)));
%   figure;plot(-xq,zq,'r-');hold on;
%   plot(-xq(idx),a.*xq(idx)+b,'k-');
end

function [a,b]= linfit(x,y)

a=0;
b=0;

if(length(x) > 0 & length(y) > 0)

X = [ones(length(x),1) x];
ab = X\y;
b=ab(1);
a=ab(2);

end
end