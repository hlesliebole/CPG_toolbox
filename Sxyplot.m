
figure;plot(wavetime(idx),Sxy(idx)+0.0025,'k.-','linewidth',2);set(gca,'fontsize',18,'linewidth',2);ylabel('Mop Sxy (m^{2})'),grid on;
title('OC478 Sxy')
set(gca,'xtick',[datetime(2006,9,11:30) datetime(2006,10,1:19)])
set(gca,'xlim',[datetime(2006,9,11) datetime(2006,10,19)])