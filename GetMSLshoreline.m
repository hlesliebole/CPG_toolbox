function [Sdatetime,xMSL]=GetMSLshoreline(MopNum)

load([ 'M'  num2str( MopNum , '%5.5i' )  'SA.mat' ],'SA');

[X1D,Z1D]=GetAllNearestPointsProfiles(SA,25,5);

xMSL=[];
xl=[X1D(1) X1D(end)];
for n=1:size(Z1D,1)
        xint=intersections(xl,[0.774 0.774],X1D,Z1D(n,:));
        if ~isempty(xint)
            xMSL(n)=max(xint);
        else
            xMSL(n)=NaN;
        end
end
Sdatetime=datetime([SA.Datenum],'convertfrom','datenum');

end