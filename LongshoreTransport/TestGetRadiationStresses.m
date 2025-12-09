
MopID='D0888';

[wavetime,Sxx,Sxy]=GetCpgMopRadiationStresses(MopID,.04,.4,0);
figure;
subplot(2,1,1)
plot(wavetime,Sxx)
hold on;
subplot(2,1,2)
plot(wavetime,Sxy)
hold on;
mean(Sxx,'omitnan')
mean(Sxy,'omitnan')

%[wavetime2,TotSxy]=GetMopSxy('D0654');

[wavetime,Sxx,Sxy]=GetCpgMopRadiationStresses(MopID,.04,.09,0);
subplot(2,1,1)
plot(wavetime,Sxx,'r-')
subplot(2,1,2)
plot(wavetime,Sxy,'r-')
mean(Sxx,'omitnan')
mean(Sxy,'omitnan')

%%
[wavetime,Sxx,Sxy]=GetCpgMopRadiationStresses(MopID,.11,.4,0);
subplot(2,1,1)
plot(wavetime,Sxx,'g-')
subplot(2,1,2)
plot(wavetime,Sxy,'g-')
mean(Sxx,'omitnan')
mean(Sxy,'omitnan')

