close all

[TZdatetime,TZ1Dmon,TZ1Dqtr,TZ1Dsea,TZ1Dann,TZ1Dglo]=...
    TimeAveragedMeanProfiles(Zdatetime,Z1Dmin);

figure;plot(X1Dcpg,TZ1Dglo,'.-')
hold on;grid on;
plot(X1Dcpg,TZ1Dmon(end,:),'r.-')
plot(X1Dcpg,TZ1Dqtr(end-1,:),'g.-')
plot(X1Dcpg,TZ1Dsea(end,:),'c.-')
plot(X1Dcpg,TZ1Dann(end,:),'k.-')

figure;hold on;grid on
for n=1:numel(TZdatetime.mon)
    plot3(TZdatetime.mon(n)+X1Dcpg*0,X1Dcpg,TZ1Dmon(n,:),'g-','linewidth',2)
end

%figure;hold on;grid on
for n=1:numel(TZdatetime.qtr)
    plot3(TZdatetime.qtr(n)+X1Dcpg*0,X1Dcpg,TZ1Dqtr(n,:),'r-','linewidth',2)
end


%figure;hold on;grid on
for n=1:numel(TZdatetime.sea)
    plot3(TZdatetime.sea(n)+X1Dcpg*0,X1Dcpg,TZ1Dsea(n,:),'c-')
end


%figure;hold on;grid on
for n=1:numel(TZdatetime.ann)
    plot3(TZdatetime.ann(n)+X1Dcpg*0,X1Dcpg,TZ1Dann(n,:),'k-','linewidth',2)
end

