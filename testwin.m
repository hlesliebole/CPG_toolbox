swin=5;
zgdsig=zgd*NaN;
for n=1:size(zgd,1) % loop through grid points with data
    zt=zgd(n,:);
    i=find(~isnan(zt));% surveys with data at this grid point
    if numel(i) >= swin % use desired running survey window if possible 
        win=swin;
    else
        win=numel(i); % otherwise use what there is
    end
    ztg=zt(i);ztsig=ztg*0+Inf; % reduce to just valid survey data time series
     for j=1:numel(i)-win+1 % loop through time series with window
         jw=j:j+win-1;
         % get sigma value for each data point
         sd=std(ztg(jw));mn=mean(ztg(jw));sig=abs((ztg(jw)-mn)./sd);
         % retain minimum sigma derived for each data point
         ztsig(jw)=min([ztsig(jw)' sig']');
     end   
     zgd(n,i)=ztsig; % save minimum sigmas for this grid point
end