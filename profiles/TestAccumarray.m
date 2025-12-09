clearvars
load M00005SA.mat
load M00005SG.mat

% find SA survey with the most gridded data
samax=0;nmax=0;
for n=1:size(SA,2)
    if numel(SA(n).Z) > nmax
        nmax=numel(SA(n).Z);
        samax=n;
    end
end

% find SG survey with the most gridded data
smax=0;nmax=0;
for n=1:size(SG,2)
    if numel(SG(n).Z) > nmax
        nmax=numel(SG(n).Z);
        smax=n;
    end
end

fprintf('Most complete grid:\n%s\n',SG(smax).File)

n=smax;
xutm=SG(n).X;yutm=SG(n).Y;z=SG(n).Z;

ColorScatter(xutm,yutm,z)


n=smax;
xutm=SA(n).X;yutm=SA(n).Y;z=SA(n).Z;
xutm=vertcat(SA.X);yutm=vertcat(SA.Y);z=vertcat(SA.Z);


Res=1; % 1m spatial resolution
% round survey x,y to desired resolution
xr=Res*round(xutm/Res); % round to Res meter spatial resolution
yr=Res*round(yutm/Res); % 

% bin and average rounded survey data by placing in unique
%  x,y data array
[ux, ~, xidx] = unique(xr);
[uy, ~, yidx] = unique(yr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray([xidx(:), yidx(:)], 1);  
%array of average of z that fall into each unique x/y combination
zavg = accumarray([xidx(:), yidx(:)], z.')./zcount;
%cmode = accumarray([xidx(:), yidx(:)], c.',[], @mode); % most common class 
%tavg = accumarray([xidx(:), yidx(:)], t.')./zcount;
%create a list of the z that fall into each unique x/y combination
%zs = accumarray([xidx(:), yidx(:)], z.', [], @(V) {V}, {});
