function d=DistanceFromLine(x,y,a,b)

% perpendicular distance of the point x,y from the line defined
%   by y = ax + b

d=abs(a.*x - y + b)./sqrt(a.^2+1);

end