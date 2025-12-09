

for mop=635:665
     
    matfile=['M' num2str(mop,'%5.5i') 'SA.mat']

load(matfile,'SA');

ndx=find( (strcmp({SA.Source},'Trk') | strcmp({SA.Source},'AtvMR')) & [SA.Datenum] < datenum(2024,1,1) );

lnorm=[];
n=0;
for m=ndx
    n=n+1;
    idx=find(SA(m).Z > 1.2 & SA(m).Z < 1.4);
    [a,b]=linfit([SA(m).X(idx)] , [SA(m).Y(idx)]);
    lnorm(n)=180-atand(a);
end

MopLnorm(mop-634)=mean(rmoutliers(lnorm),'omitnan');
MopLnormSdev(mop-634)=std(rmoutliers(lnorm),'omitnan');

end

save SolanaLidarShoreNormals.mat MopLnorm MopLnormSdev
figure;errorbar(635:665,MopLnorm,MopLnormSdev)

