function [MinRmsd,xlag,xlength]=GetRmsdGMP(GMP,mop1,mop2,minZ,maxZ)

MinRmsd=[];
xlag=[];
xlength=[];

imop1=find([GMP.Mop] == mop1);
imop2=find([GMP.Mop] == mop2);
MinRmsd=Inf;
for dx=-50:50
xp1=round(GMP(imop1).X);
xp2=round(GMP(imop2).X+dx);
gp1=GMP(imop1).Zglobal;
gp2=GMP(imop2).Zglobal;
gp1(gp1 < minZ)=NaN;
gp1(gp1 > maxZ)=NaN;
gp2(gp2 < minZ)=NaN;
gp2(gp2 > maxZ)=NaN;
idx1=find(~isnan(gp1));
idx2=find(~isnan(gp2));

if ~isempty(idx1) && ~isempty(idx2)
    
overlap=max([xp1(idx1(1)) xp2(idx2(1))]):min([xp1(idx1(end)) xp2(idx2(end))]);

gp1=gp1(ismember(xp1,overlap));
gp2=gp2(ismember(xp2,overlap)); 

diff=gp2-gp1;
rmsd=sqrt(sum(diff.^2,'omitnan')/numel(~isnan(diff)));
%fprintf('%6.1f  %8.6f  %i\n',dx,rmsd,numel(~isnan(diff)))
MinRmsd=min([MinRmsd rmsd]);
if rmsd == MinRmsd;xlag = dx;xlength=numel(overlap);end
end

end
end