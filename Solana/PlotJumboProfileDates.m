load M00653SA.mat

figure('position',[ 140         313        1096         372]);

idx=find(contains({SA.File},'umbo') | contains({SA.File},'etski') );
ng=0;
for n=idx
    if min(SA(n).Z) < -3
        hold on;
        ng=ng+1;
      p2=plot([SA(n).Datenum],[SA(n).Datenum]*0+653,'ro','linewidth',2);
    end
end

load M00654SA.mat
idx=find(contains({SA.File},'umbo') | contains({SA.File},'etski') );


ng1=0;
for n=idx
    if min(SA(n).Z) < -3
        hold on;
         ng1=ng1+1;
      p1=plot([SA(n).Datenum],[SA(n).Datenum]*0+654,'b+','linewidth',2);
    end
end


grid on;box on;
datetick('x')
legend([p1 p2],['Mop 654; N=' num2str(ng1)],['Mop 653; N=' num2str(ng)],'fontsize',18)
set(gca,'ylim',[650 658],'ytick',[653 654],'fontsize',16);

title('Jumbo Survey History (blue plus only = 200m spacing even Mops survey; w/ red circle = 100m spacing')
makepng('SolanaJumboSurveyHistory.png')