load TP5mTotalFlux.mat;
figure('position',[454   427   825   336]);
FhwTS(FhwTS < .28)=0;
esum=cumsum(-(FhwTS(:,501).^3));C=-23/esum(end);esum=esum*C;
plot(datenum(xt),esum);datetick