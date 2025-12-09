MopNumber=582;
load VosTransects.mat
for n=1:size(VosT,2)
    if isempty(VosT(n).BackMop)
       VosT(n).BackMop=NaN;
       VosT(n).BackMopXoffset=NaN;
    end
end

    idx=find(round([VosT.BackMop]) == MopNumber); 
    

if ~isempty(idx)

 if numel(idx) > 1
  fprintf('More than one Vos Transect is closest to Mop: %i\n',MopNumber)
   for i=1:numel(idx)
       VosT(idx(i))
   end
 end
 
  url=['http://coastsat.wrl.unsw.edu.au/time-series/' VosT(idx(1)).Name '/'];
  s=webread(url);

else
    fprintf('No Vos Transect is closest to Mop: %i\n',MopNumber)
end

ss=strsplit(s,'\n');
k=0;
for n=1:size(ss,2)
    if ~isempty(ss{n})
       st=strsplit(regexprep(ss{n},',',' '),' ');
       k=k+1;
       dt(k)=datetime(st{1});
       xMSL(k)=str2double(st{3});
       xMSLmop(k)=str2double(st{3})+VosT(idx(1)).BackMopXoffset;
    end   
end

