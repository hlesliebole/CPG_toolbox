function [a,b]= linfit(x,y)

a=0;
b=0;

if(length(x) > 0 & length(y) > 0)

X = [ones(length(x),1) x];
ab = X\y;
b=ab(1);
a=ab(2);

end