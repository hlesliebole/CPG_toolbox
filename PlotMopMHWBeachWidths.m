
Mopnum=848;
[Sdatetime,xMHW]=GetMHWshoreline(Mopnum);
figure('position',[  33  380  1275   364]);
plot(Sdatetime,xMHW,'*-');ylabel('MHW Beach Width');...
    title(['Mop ' num2str(Mopnum)]);set(gca,'fontsize',18)
grid on
set(gca,'xtick',datetime(2000:1:2025,1,1))