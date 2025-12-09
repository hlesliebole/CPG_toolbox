
dz=0.1;
zmax=dz*ceil(max(Z(:))/dz);
zmin=dz*floor(min(Z(:))/dz);

[nz,zn]=histcounts(Z(:),zmin:dz:zmax);
zn=zn(1)+dz/2:dz:zn(end-1)+dz/2;

figure;
plot(zn,nz,'k-',zn,nz,'k.','markersize',10);grid on
xlabel('Bin Elevation (m, NAVD88)');
ylabel('Number of Grid Points in Bin');
hold on;

[dmin,dmax]=GetQuantCutoffs(Z,1);

plot([dmin dmin],[0 max(nz)],'k--');
plot([dmax dmax],[0 max(nz)],'k--');


[sx,sy]=gradient(Z);
% use beta = max gradient component at each x,y location
beta=max(cat(3,abs(sx),abs(sy)),[],3);

% loop through elevation levels and find the one
%   with at least 200 points that has the best
%   (smallest standard deviation) linear fit line.

minstd=Inf;abest=[];bbest=[];zbest=[];
for z=zmin:dz:zmax-dz
    idx=find(Z(:) >= z & Z(:) < z+dz); % points in elev range
    if length(idx) > 200 % skip if too few elevs in this range
        [a,b]=linfit(X(idx),Y(idx)); % best fit line to elv points
        d=(a.*X(idx) - Y(idx) + b)./sqrt(a.^2+1); % distance from line
        if minstd > nanstd(d) % standard dev sigma distance
            minstd=nanstd(d);abest=a;bbest=b;zbest=z+dz/2;
        end
    end
end

% best fit normal (True compass coming from) 
BchNorm=360-atand(abest);

yyaxis right
smed=[];
for z=zmin:dz:zmax-dz
    idx=find(Z(:) >= z & Z(:) < z+dz);
    if ~isempty(idx) % skip if no elevs in this range
    plot((z+dz/2)*ones(length(idx),1),s(idx),'g.');
    plot((z+dz/2),nanmedian(s(idx)),'b.');
    smed=[smed nanmedian(s(idx))];
    [dmin,dmax]=GetQuantCutoffs(s(idx),2);
    plot((z+dz/2),dmin,'r.');
    plot((z+dz/2),dmax,'r.');
    end
end
ylabel('Grid Point Gradients');
%figure;plot(zn,smed)
