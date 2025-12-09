x = [14.9 1.7 0.0 10.9 0.0];
y = [11.3 9.1 23.7 12.8 2.9];
z = [5.32787E-17 2.93234E-16 2.09997E-16 5.45E-17 4.55E-16];
B = [x(:) y(:) ones(size(x(:)))] \ z(:);
xv = linspace(min(x), max(x), 10)';
yv = linspace(min(y), max(y), 10)';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);
scatter3(x,y,z, 'filled')
hold on
mesh(X, Y, Z, 'FaceAlpha', 0.5)
hold off
view(-120, 35)
title(sprintf('Z = %+.3E\\cdotX %+.3E\\cdotY %+.3E', B)) 
atan2d(max(diff(Y(:,1))),max(diff(X(1,:))))