% Mop number to grab CoastSat data for
MopNumber=582;

% load struct array containing info about the relationship of each
%  Mop transect to a CoastSat transect(s).
load VosTransects2.mat

% set null values in struct array to NaNs
for n=1:size(VosT,2)
    if isempty(VosT(n).BackMop)
       VosT(n).BackMop=NaN;
       VosT(n).BackMopXoffset=NaN;
    end
end

% find the CoastSat transect closest to the Mop transect
idx=find(round([VosT.BackMop]) == MopNumber); 
    

if ~isempty(idx)

 if numel(idx) > 1
  fprintf('More than one Vos Transect is closest to Mop: %i\n',MopNumber)
   for i=1:numel(idx)
       VosT(idx(i))
   end
 end
 
 % get CoastSat transect time series of shoreline location
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
       xShoreline(k)=str2double(st{3});
       % transform to a Mop transect location
       xMopShoreline(k)=str2double(st{3})+VosT(idx(1)).BackMopXoffset;
    end   
end

% Make example plot
figure('position',[ 98         190        1146         540]);
plot(dt,xShoreline,'bo-','DisplayName',...
    ['CoastSat ' VosT(idx(1)).Name ' Transect Shoreline Distance (m) from Back Point']);
hold on;grid on;
plot(dt,xMopShoreline,'r*-','DisplayName',...
    ['Mop Transformed CoastSat Shoreline Distance (m) from Mop ' num2str(MopNumber) ' Back Point']);
xlabel('Date');ylabel('Cross-shore Distance (m)');set(gca,'fontsize',14)
legend('location','northoutside','interpreter','none')
