%function [x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,NumTrans,ExtendLine)

load MopTableUTM.mat

%  Used by BuildSMmatfiles.m

%  Calculates the interpolated 1m xshore resolution transect line, and NumTrans 
%   subtransect lines, for/within the MopNumber survey area bounds. Extends
%   the lines -/+ ExtendLine(1)/ ExtendLine(2) meters from the 
%   back beach (-) and offshore (+) mop line, using the Mop table for
%   reference.

% eg. for 20 subtransects that extend -100m back from back beach
%      mop line:
%
%  [x1d,xt,yt,xst,yst]=GetTransectLines(Mop,580,20,[-100 0])

%  returns transects all of the same length = longest of the
%  subtransects


%-------------------------------------------------------------
% first, the subtransect lines
%-------------------------------------------------------------

 % tweak offshore mop point locations to be at .999 of the distance of the
 % actual mop transect to prevent neighboring mops sharing the exact same 
 % offshore mop location.
 
 for m=MopNumber-1:MopNumber+1
  xa=Mop.BackXutm(m)+0.999*(Mop.OffXutm(m)-Mop.BackXutm(m));
  xyratio=(xa-Mop.BackXutm(m))/(Mop.OffXutm(m)-Mop.BackXutm(m));
  ya=Mop.BackYutm(m)+xyratio*(Mop.OffYutm(m)-Mop.BackYutm(m));
  Mop.OffXutm(m)=xa;
  Mop.OffYutm(m)=ya;
 end
 
 % back beach neighboring mop midpoints
 xb1=mean([Mop.BackXutm(MopNumber-1) Mop.BackXutm(MopNumber)]);
 yb1=mean([Mop.BackYutm(MopNumber-1) Mop.BackYutm(MopNumber)]);
 xb2=Mop.BackXutm(MopNumber);
 yb2=Mop.BackYutm(MopNumber);
 xb3=mean([Mop.BackXutm(MopNumber) Mop.BackXutm(MopNumber+1)]);
 yb3=mean([Mop.BackYutm(MopNumber) Mop.BackYutm(MopNumber+1)]);
 % offshore neighboring mop midpoints
 xo1=mean([Mop.OffXutm(MopNumber-1) Mop.OffXutm(MopNumber)]);
 yo1=mean([Mop.OffYutm(MopNumber-1) Mop.OffYutm(MopNumber)]);
 xo2=Mop.OffXutm(MopNumber);
 yo2=Mop.OffYutm(MopNumber);
 xo3=mean([Mop.OffXutm(MopNumber) Mop.OffXutm(MopNumber+1)]);
 yo3=mean([Mop.OffYutm(MopNumber) Mop.OffYutm(MopNumber+1)]);
 
 nst=NumTrans; % number of subtransects 
 % increase nst by 2 when calculating x,y step size to keep
 %  first and last subtransects inside mop area boundaries
 %  for 2d interpolation. 
%  xbst=linspace(xb1,xb2,nst+2);xbst=xbst(2:end-1);
%  ybst=linspace(yb1,yb2,nst+2);ybst=ybst(2:end-1);
%  xost=linspace(xo1,xo2,nst+2);xost=xost(2:end-1);
%  yost=linspace(yo1,yo2,nst+2);yost=yost(2:end-1);
 
 [xbst,ybst]=EqualSpacedPoints([xb1 xb2 xb3],[yb1 yb2 yb3],NumTrans);
 xo1 
 xo2 
 xo3
 yo1
 yo2
 yo3
 
 [xost,yost]=EqualSpacedPoints([xo1 xo2 xo3],[yo1 yo2 yo3],NumTrans);
 
%  270-atan2d(yb2-yo2,xb2-xo2)
%  270-atan2d(yb2-yo2+2,xb2-xo2)
%  270-atan2d(yb2-yo2-2,xb2-xo2)
%  270-atan2d(ybst-yost,xbst-xost)
%  mean(270-atan2d(ybst-yost,xbst-xost))
 
 % get max subtransect length
 dmax=0;
 for n=1:nst
     dist=pdist([xbst(n),ybst(n);xost(n),yost(n)]);
     dmax=max([dmax dist]);
 end
 
 % include the actual Mop transect distance 
 x1=Mop.BackXutm(MopNumber);y1=Mop.BackYutm(MopNumber);
 x2=Mop.OffXutm(MopNumber);y2=Mop.OffYutm(MopNumber);
 dist=pdist([x1,y1;x2,y2]);
 
 dmax=max([dmax dist]); % max of the transect & subtransects lengths
 dmax=ceil(dmax); % round up to whole meter
 

 % now loop through subtransect endpoints and make 1m xshore interpolated
 %  transects
 for n=1:nst
   x1=xbst(n);y1=ybst(n);x2=xost(n);y2=yost(n);
   ang=atan2(y2-y1,x2-x1); % radian angle between back and offshore point
   dist=pdist([x1,y1;x2,y2]); % distance between back and offshore point
   % adjust offshore extension so all transects have same xshore length
   ExtendOff=ExtendLine(2)+(dmax-dist);
 
   % 1m spaced points along transect 
   % x coords of 1m spaced points along transect
   xst(n,:)=x1+ExtendLine(1)*cos(ang):cos(ang):x2+ExtendOff*cos(ang);
   % y coords of 1m spaced points along transect
   yst(n,:)=y1+ExtendLine(1)*sin(ang):sin(ang):y2+ExtendOff*sin(ang);  
 
 end
 
 % xshore distance m from back beach point
 x1d=ExtendLine(1):dmax+ExtendLine(2);
 
%-------------------------------------------------------------
% now the single Mop transect line
%-------------------------------------------------------------
    
 % Mop transect 1m spaced xshore points
 x1=Mop.BackXutm(MopNumber);y1=Mop.BackYutm(MopNumber);
 x2=Mop.OffXutm(MopNumber);y2=Mop.OffYutm(MopNumber);
 ang=atan2(y2-y1,x2-x1); % radian angle between back and offshore point
 dist=pdist([x1,y1;x2,y2]); % distance between back and offshore point
 % adjust offshore extension so all transects have same xshore length
 ExtendOff=ExtendLine(2)+(dmax-dist);
 
 % x coords of 1m spaced points along transect
 xt=x1+ExtendLine(1)*cos(ang):cos(ang):x2+ExtendOff*cos(ang);
 % y coords of 1m spaced points along transect
 yt=y1+ExtendLine(1)*sin(ang):sin(ang):y2+ExtendOff*sin(ang);  
 
end