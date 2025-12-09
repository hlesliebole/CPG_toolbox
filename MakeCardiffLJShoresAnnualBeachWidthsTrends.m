% This version plots linear best fit trend line of 
% beach widths between el nino years.

% starts with 2004, requires at least 3 values to
% show a trend line, otherwise a dotted line
% between two values, and no line for 0-1 values. 

% It makes two plots, one with beaches where the most recent trend line
%  (post 2016 el nino) is positive, and one with beaches that have 
%  no data/line or a negative trend.
addpath '/Users/William/Desktop/ShoreBox/v2.0'
addpath '/Users/William/Desktop/ShoreBox/v1.0'

% ReportYear=2019;
% SurveyFraction=0.7;

MakeStateBeachMapFigure2
set(gca,'position',[0.01 0.005 0.33 0.95]); % adjust map axes position
set(gca,'xlim',[-117.3787 -117.1890],'ylim',[32.8069   33.0361]);
delete(sbb(5:end));delete(sbl(5:end));delete(sbt(5:end))
sbt(3).Position=[-117.3097+0.025 32.9108 0];
sbl(3).XData=[-117.2597 -117.3097+0.025];
sbt(3).FontSize=18;
sbt(4).Position=[-117.3316+0.025 33.0054 0];
sbl(4).XData=[-117.2816 -117.3316+0.025];
sbt(4).FontSize=18;
xy =[
 -117.2840   32.9959
 -117.2762   32.9976
 -117.2704   32.9770
 -117.2778   32.9754
 -117.2836   32.9959];
p=plot(xy(:,1),xy(:,2),'g-','linewidth',2);
t1=text(-117.3847,32.985,'Solana Beach','fontsize',18,'fontweight','bold','color',[.95 .95 .95]);
pl=plot([-117.3037 -117.2754],[32.9844 32.9844],'-','color',[.95 .95 .95]);
    
xy =[
 -117.2766   32.9717
 -117.2688   32.9737
 -117.2643   32.9520
 -117.2721   32.9506
 -117.2766   32.9721];
p=plot(xy(:,1),xy(:,2),'g-','linewidth',2);
t2=text( -117.3365,32.9585,'Del Mar','fontsize',18,'fontweight','bold','color',[.95 .95 .95]);
pl2=plot([-117.290 -117.2700],[32.9610   32.9610],'-','color',[.95 .95 .95]);

xy =[
 -117.2598   32.8687
 -117.2531   32.8654
 -117.2570   32.8523
 -117.2638   32.8531
 -117.2598   32.8687];
p=plot(xy(:,1),xy(:,2),'g-','linewidth',2);
t3=text(-117.3647,32.8593,'La Jolla Shores','fontsize',18,'fontweight','bold','color',[.95 .95 .95]);
pl2=plot([-117.2732 -117.2605],[32.8593   32.8593],'-','color',[.95 .95 .95]);




set(gcf,'position',[15         343        1365         430])


%A2DefineStateBeaches

%lt=[':o';':^';':x';':d';':v';':*';':p';':s';':h';':>';':+';':<'];
lt=['o';'^';'x';'d';'v';'*';'p';'s';'h';'>';'+';'<'];

%figure('position',[440   189   800   600]);
ax2=axes('position',[.37 .15 .58 .65],'yticklabels',[]);hold on;


set(gca,'ylim',[10 70]);
yl=get(gca,'ylim');
% highlight el nino years
ey=1998;

fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
ey=2003;
fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
ey=2010;
fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
ey=2016;
fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);

% major nourishments
ms=2001+4/12;me=2001+9/12; % regional beach sand project I
mn(1)=fill([ms me  me ms ms],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.9 0.9 0.],'facealpha',0.3,'edgecolor',[0.9 0.9 0.0]);
ms=2012+9/12;me=2012+11/12; % regional beach sand project II
mn(2)=fill([ms me  me ms ms],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0. 0.9 0.],'facealpha',0.3,'edgecolor',[0. 0.9 0.]);
ms=2018+2/12;me=2018+6/12; % san elijo lagoon restoration project
mn(3)=fill([ms me  me ms ms],...
    [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.0 0.9 0.9],'facealpha',0.3,'edgecolor',[0.0 0.9 0.9]);

 nn=1;
 dname='Cardiff';
[xmbw,mbw]=GetReachAnnualMeanBeachWidths([666 682]);

    h(nn)=plot(xmbw+.3,mbw,'bx','linewidth',2,'markersize',10,'DisplayName',dname);hold on;
    
    j=find(xmbw > 2015 & xmbw < 2021);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        col=get(h(nn),'color'); 
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end
    j=find(xmbw > 1997 & xmbw < 2003);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  
    j=find(xmbw > 2002 & xmbw < 2010);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  
    j=find(xmbw > 2009 & xmbw < 2016);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  

    nn=2;
 dname='Solana';
[xmbw,mbw]=GetReachAnnualMeanBeachWidths([639 663]);

    h(nn)=plot(xmbw+.3,mbw,'r*','linewidth',2,'markersize',10,'DisplayName',dname);hold on;
    
    j=find(xmbw > 2015 & xmbw < 2021);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        col=get(h(nn),'color'); 
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end
    j=find(xmbw > 1997 & xmbw < 2003);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  
    j=find(xmbw > 2002 & xmbw < 2010);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  
    j=find(xmbw > 2009 & xmbw < 2016);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  

    nn=3;
 dname='North Del Mar';
[xmbw,mbw]=GetReachAnnualMeanBeachWidths([621 635]);

    h(nn)=plot(xmbw+.3,mbw,'gs','linewidth',2,'markersize',10,'DisplayName',dname);hold on;
    set(h(nn),'color',[0 .7 0]);
    j=find(xmbw > 2015 & xmbw < 2021);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        col=get(h(nn),'color'); 
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end
    j=find(xmbw > 1997 & xmbw < 2003);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  
    j=find(xmbw > 2002 & xmbw < 2010);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  
    j=find(xmbw > 2009 & xmbw < 2016);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  
    
    nn=4;
 dname='South Torrey';
[xmbw,mbw]=GetReachAnnualMeanBeachWidths([540 565]);

    h(nn)=plot(xmbw+.3,mbw,'m^','linewidth',2,'markersize',10,'DisplayName',dname);hold on;
    
    j=find(xmbw > 2015 & xmbw < 2021);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        col=get(h(nn),'color'); 
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end
    j=find(xmbw > 1997 & xmbw < 2003);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  
    j=find(xmbw > 2002 & xmbw < 2010);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  
    j=find(xmbw > 2009 & xmbw < 2016);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end 
    
    nn=5;
 dname='La Jolla Shores';
[xmbw,mbw]=GetReachAnnualMeanBeachWidths([497 514]);

    h(nn)=plot(xmbw+.3,mbw,'ko','linewidth',2,'markersize',10,'DisplayName',dname);hold on;
    
    j=find(xmbw > 2015 & xmbw < 2021);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        col=get(h(nn),'color'); 
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end
    j=find(xmbw > 1997 & xmbw < 2003);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  
    j=find(xmbw > 2002 & xmbw < 2010);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  
    j=find(xmbw > 2009 & xmbw < 2016);
    if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
        xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
    end  
    
set(gca,'fontsize',16)
xlabel('Beach Year');
box on

set(gca,'xtick',[1994:2:str2num(datestr(date,'yyyy'))+1])
%xlabel('Beach Year');
set(gca,'ylim',[10 70]);
yl=get(gca,'ylim');
yyaxis right;
set(gca,'ylim',yl,'ycolor','k');
ylabel('Mean Beach Width (m)');
set(ax2,'xlim',[1998 2021]);
grid on
legend(h(1:5),'location','west')
title('Annual Mean MHW Beach Width Change','fontsize',22)
ax3=axes('position',[.392 .77 .1 .1]);
l3=legend(ax3,mn,'Regional Beach Sand I','Regional Beach Sand II','San Elijo Lagoon Restoration','location','northeast',...
    'fontsize',12,'linewidth',2,'edgecolor','k');title(l3,'Major Nourishment Projects')
axis off

ax4=axes('position',[.824 .805 .1 .1]);
plot([1998],[0],'-','linewidth',10,'color',[.8 .8 .8],'DisplayName','El Nino Years');
legend('fontsize',13,'location','south','linewidth',2,'fontweight','bold');
axis off

makepng('CardiffToLJshoresMhwBeachWidthChange.png')


% ax1=axes('position',[.37 .55 .58 .4],'yticklabels',[]);hold on;
% set(gca,'ylim',[10 80]);
% yl=get(gca,'ylim');
% % highlight el nino years
% ey=1998;
% 
% fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
%     [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
% ey=2003;
% fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
%     [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
% ey=2010;
% fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
%     [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
% ey=2016;
% fill([ey-.25 ey+.75 ey+.75 ey-.25 ey-.25],...
%     [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.8 0.8 0.8],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8]);
% 
% % major nourishments
% ms=2001+4/12;me=2001+9/12; % regional beach sand project I
% mn(1)=fill([ms me  me ms ms],...
%     [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.9 0.9 0.],'facealpha',0.3,'edgecolor',[0.9 0.9 0.0]);
% ms=2012+9/12;me=2012+11/12; % regional beach sand project II
% mn(2)=fill([ms me  me ms ms],...
%     [yl(1) yl(1) yl(2) yl(2) yl(1)],[0. 0.9 0.],'facealpha',0.3,'edgecolor',[0. 0.9 0.]);
% ms=2018+2/12;me=2018+6/12; % san elijo lagoon restoration project
% mn(3)=fill([ms me  me ms ms],...
%     [yl(1) yl(1) yl(2) yl(2) yl(1)],[0.0 0.9 0.9],'facealpha',0.3,'edgecolor',[0.0 0.9 0.9]);
% 


%ax2=axes('position',[.37 .1 .58 .84]);hold on;
%set(gca,'color','none');
% 
% nn=0;
% for n=size(StateBeach,2):-1:1
%     
% if ~isempty(StateBeach(n).South)
%     
%     nn=nn+1;
%     [xmbw,mbw]=getAnnualMeanBeachWidths(StateBeach(n).name,StateBeach(n).North,SurveyFraction);
%     %dname=['North ' StateBeach(n).title]; 
%     dname=strrep(['North ' StateBeach(n).title],' State',''); 
%     dname=strrep(dname,' Beach','');
%     dname=strrep(dname,' Park','');
%     
%     j=find(xmbw > 2015);
%     if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
%         % force border south to be accreting axis because there is no data
%         if(a > 0 || n == 1);axes(ax2);ax(nn)=2;else;axes(ax1);ax(nn)=1;end
%         
%     h(nn)=plot(xmbw+.3,mbw,lt(nn,:),'linewidth',2,'markersize',10,'DisplayName',dname);hold on;
%     if(nn == 6);set(h(nn),'color','m');end
%     col=get(h(nn),'color');%set(h(nn),'markerfacecolor',col);     
%         xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
%     
%     
%     else
%         % force border south to be accreting axis because there is no data
%         if(n == 1);axes(ax2);ax(nn)=2;else;axes(ax1);ax(nn)=1;end
%         h(nn)=plot(xmbw+.3,mbw,lt(nn,:),'linewidth',2,'markersize',10,'DisplayName',dname);hold on;
%         if(nn == 6);set(h(nn),'color','m');end
%         col=get(h(nn),'color');%set(h(nn),'markerfacecolor',col);   
%     end 
%     j=find(xmbw > 1997 & xmbw < 2003);
%     if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
%         xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
%     end  
%     j=find(xmbw > 2002 & xmbw < 2010);
%     if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
%         xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
%     end  
%     j=find(xmbw > 2009 & xmbw < 2016);
%     if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
%         xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
%     end  
%     
%      
%     nn=nn+1;
%     [xmbw,mbw]=getAnnualMeanBeachWidths(StateBeach(n).name,StateBeach(n).South,SurveyFraction);
%     %dname=['South ' StateBeach(n).title]; 
%     dname=strrep(['South ' StateBeach(n).title],' State',''); 
%     dname=strrep(dname,' Beach','');
%     dname=strrep(dname,' Park','');
%      
%     j=find(xmbw > 2015);
%     if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
%         if(a > 0);axes(ax2);ax(nn)=2;else;axes(ax1);ax(nn)=1;end
%         
%     h(nn)=plot(xmbw+.3,mbw,lt(nn,:),'linewidth',2,'markersize',10,'DisplayName',dname);hold on;
%     if(nn == 6);set(h(nn),'color','m');end
%     col=get(h(nn),'color');%set(h(nn),'markerfacecolor',col);     
%         xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
%     
%     elseif (length(j) == 1 || n ==1 ) % Border south is special case
%         axes(ax2);ax(nn)=2;
%         h(nn)=plot(xmbw+.3,mbw,lt(nn,:),'linewidth',2,'markersize',10,'DisplayName',dname);hold on;
%         if(nn == 6);set(h(nn),'color','m');end
%         col=get(h(nn),'color');%set(h(nn),'markerfacecolor',col); 
%     else
%         axes(ax1);ax(nn)=1;
%         h(nn)=plot(xmbw+.3,mbw,lt(nn,:),'linewidth',2,'markersize',10,'DisplayName',dname);hold on;
%         if(nn == 6);set(h(nn),'color','m');end
%         col=get(h(nn),'color');%set(h(nn),'markerfacecolor',col);   
%     end 
%     j=find(xmbw > 1997 & xmbw < 2003);
%     if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
%         xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
%     end  
%     j=find(xmbw > 2002 & xmbw < 2010);
%     if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
%         xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
%     end  
%     j=find(xmbw > 2009 & xmbw < 2016);
%     if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
%         xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
%     end  
%     
%     
% else
%     
%     nn=nn+1;
%     [xmbw,mbw]=getAnnualMeanBeachWidths(StateBeach(n).name,[],SurveyFraction);
%     %dname=[StateBeach(n).title];
%     dname=strrep(StateBeach(n).title,' State',''); 
%     dname=strrep(dname,' Beach','');
%     dname=strrep(dname,' Park','');
%     
%     
%     j=find(xmbw > 2015);
%     if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
%         if(a > 0);axes(ax2);ax(nn)=2;else;axes(ax1);ax(nn)=1;end
%         
%     h(nn)=plot(xmbw+.3,mbw,lt(nn,:),'linewidth',2,'markersize',10,'DisplayName',dname);hold on;
%     if(nn == 6);set(h(nn),'color','m');end
%     col=get(h(nn),'color');%set(h(nn),'markerfacecolor',col);     
%         xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
%     else
%         axes(ax1);ax(nn)=1;
%         h(nn)=plot(xmbw+.3,mbw,lt(nn,:),'linewidth',2,'markersize',10,'DisplayName',dname);hold on;
%         if(nn == 6);set(h(nn),'color','m');end
%         col=get(h(nn),'color');%set(h(nn),'markerfacecolor',col);   
%     end 
%     j=find(xmbw > 1997 & xmbw < 2003);
%     if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
%         xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
%     end  
%     j=find(xmbw > 2002 & xmbw < 2010);
%     if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
%         xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
%     end  
%     j=find(xmbw > 2009 & xmbw < 2016);
%     if(length(j) > 1);[a,b]=linfit(xmbw(j)',mbw(j)');
%         xy=min(xmbw(j)):max(xmbw(j));plot(xy+.3,a.*xy+b,'-','color',col,'linewidth',2);
%     end  
%     
%     
% end
%     
% % matfile=[StateBeach(n).name 'AnnualBeachWidth.mat'];
% % eval(['load ' matfile ]);
% %plot(xmbw,mbw,lt(n,:),'linewidth',2,'markersize',10);hold on;
% %plot(xmbw,mbw,'.','markersize',15);
% 
% end 
% %axes(ax2)
% 
% 
% 
% axes(ax1);
% grid on
% grid minor
% l1=legend(h(ax == 1),'fontsize',14,'location','northwest','linewidth',5,'edgecolor','r');hold on; 
% title(l1,{'Eroding Beaches','\rm 2016-Present'});
% %plot([1998],[0],'-','linewidth',10,'color',[.8 .8 .8],'DisplayName','El Nino Years');
% pos=get(l1,'position');pos(1)=pos(1)-0.8*pos(3);set(l1,'position',pos);
% pos=get(l1,'position');pos(2)=pos(2)-0.05*pos(4);set(l1,'position',pos);
% %legend(StateBeach(size(StateBeach,2):-1:1).title,'fontsize',14,'location','northwest')
% set(gca,'fontsize',16)
% %title('Annual Mean Beach Widths','fontsize',18);
% title('San Diego County State Beaches: Annual Mean Beach Width Change','fontsize',18);
% 
% set(gca,'xtick',[1994:2:str2num(datestr(date,'yyyy'))+1])
% set(gca,'ylim',[10 80]);
% yl=get(gca,'ylim');
% yyaxis right;
% set(gca,'ylim',yl,'ycolor','k');
% ylabel('Mean Beach Width (m)');
% set(ax1,'xlim',[1996 2021]);
% box on
% 
% axes(ax2);
% grid on
% grid minor
% l2=legend(h(ax == 2),'fontsize',14,'location','northwest','linewidth',5,'edgecolor','g');hold on; 
% title(l2,{'Accreting Beaches','\rm 2016-Present'});
% pos=get(l2,'position');pos(1)=pos(1)-0.77*pos(3);set(l2,'position',pos);
% pos=get(l2,'position');pos(2)=pos(2)-.05*pos(4);set(l2,'position',pos);
% %plot([1998],[0],'-','linewidth',10,'color',[.8 .8 .8],'DisplayName','El Nino Years');
% %legend(StateBeach(size(StateBeach,2):-1:1).title,'fontsize',14,'location','northwest')
% set(gca,'fontsize',16)
% xlabel('Beach Year');
% box on
% 
% set(gca,'xtick',[1994:2:str2num(datestr(date,'yyyy'))+1])
% %xlabel('Beach Year');
% set(gca,'ylim',[10 80]);
% yl=get(gca,'ylim');
% yyaxis right;
% set(gca,'ylim',yl,'ycolor','k');
% ylabel('Mean Beach Width (m)');
% set(ax2,'xlim',[1996 2021]);
% 

% 
% SIOCPG=imread('SIOCPGlogo8.png');
% ax5=axes('position',[.0 .945 .35 .056]);
% image(SIOCPG);
% set(ax5,'dataaspectratio',[1 1 1]);
% axes(ax5);axis off;
% %print(gcf,'-dpng','-r50','-loose','All.png');
% 
% % %ryear=xmbw(end);
% pngfile=[num2str(ReportYear) 'AnnualReportAllBeachWidths.png']
% print(gcf,'-dpng','-r80','-loose',pngfile);
% 
% save BeachReportColors.mat h 