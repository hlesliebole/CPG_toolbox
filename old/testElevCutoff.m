

q=.90:0.001:1;
tfmax=0;fac=0;
while tfmax == 0
fac=fac+1;
tf=ischange(fac*gradient(quantile(Z(:),q)),'linear');
tfmax=max(tf);
end

topcutoff=quantile(Z(:),q(find(tf == 1,1,'first')));

figure;subplot(2,2,2);plot(q,quantile(Z(:),q));hold on;
p1=plot([q(1) q(end)],[topcutoff topcutoff],'k--');
xlabel('quantile');ylabel('z(m)');
legend(p1,'cutoff','location','northwest');
title('Extreme High Elevations');grid on;

q=0.0:0.001:.1;
tfmax=0;fac=0;
while tfmax == 0
fac=fac+1;
tf=ischange(fac*gradient(quantile(Z(:),q)),'linear');
tfmax=max(tf);
end

botcutoff=quantile(Z(:),q(find(tf == 1,1,'last')));

subplot(2,2,1);plot(q,quantile(Z(:),q));hold on;
p1=plot([q(1) q(end)],[botcutoff botcutoff],'k--');
xlabel('quantile');ylabel('z(m)');
legend(p1,'cutoff','location','southeast');
title('Extreme Low Elevations');grid on;

%-----

q=.95:0.001:1;
tfmax=0;fac=0;
while tfmax == 0
fac=fac+10;
tf=ischange(fac*(quantile(fx(:),q)),'linear');
tfmax=max(tf);
end

topcutoff=quantile(fx(:),q(find(tf == 1,1,'first')));

subplot(2,2,4);plot(q,quantile(fx(:),q));hold on;
p1=plot([q(1) q(end)],[topcutoff topcutoff],'k--');
xlabel('quantile');ylabel('slope');
legend(p1,'cutoff','location','northwest');
title('Extreme High Slopes');grid on;

q=0.0:0.001:.05;
tfmax=0;fac=0;
while tfmax == 0
fac=fac+10;
tf=ischange(fac*(quantile(fx(:),q)),'linear');
tfmax=max(tf);
end

botcutoff=quantile(fx(:),q(find(tf == 1,1,'last')));

subplot(2,2,3);plot(q,quantile(fx(:),q));hold on;
p1=plot([q(1) q(end)],[botcutoff botcutoff],'k--');
xlabel('quantile');ylabel('slope');
legend(p1,'cutoff','location','southeast');
title('Extreme Low Slopes');grid on;

%-----

q=.95:0.001:1;
tfmax=0;fac=0;
while tfmax == 0
fac=fac+10;
tf=ischange(fac*(quantile(fy(:),q)),'linear');
tfmax=max(tf);
end

topcutoff=quantile(fy(:),q(find(tf == 1,1,'first')));

subplot(2,2,4);plot(q,quantile(fy(:),q),'r-');hold on;
p1=plot([q(1) q(end)],[topcutoff topcutoff],'r--');
xlabel('quantile');ylabel('slope');
legend(p1,'Y slope cutoff','location','northwest');
title('Extreme High Slopes');grid on;

q=0.0:0.001:.05;
tfmax=0;fac=0;
while tfmax == 0
fac=fac+10;
tf=ischange(fac*(quantile(fy(:),q)),'linear');
tfmax=max(tf);
end

botcutoff=quantile(fy(:),q(find(tf == 1,1,'last')));

subplot(2,2,3);plot(q,quantile(fy(:),q),'r-');hold on;
p1=plot([q(1) q(end)],[botcutoff botcutoff],'r--');
xlabel('quantile');ylabel('slope');
legend(p1,'Y slope cutoff','location','southeast');
title('Extreme Low Slopes');grid on;

    
