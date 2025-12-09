% makes a plot showing the goodness QC sensitivity 
% the input histogram bin size for a set of x,y,z survey data.

fprintf(' Running survey QC for 26 different histogram bin sizes...\n') 
ng0=[];ng1=[];ng2=[];ng3=[];ng4=[];ng5=[];
for dz=0.005:.001:0.03
  g=ElevGoodnessQCv1(x,y,z,dz);
  ng0=[ng0 length(find(g(:) == 0))];
  ng1=[ng1 length(find(g(:) == 1))];
  ng2=[ng2 length(find(g(:) == 2))];
  ng3=[ng3 length(find(g(:) == 3))];
  ng4=[ng4 length(find(g(:) == 4))];
  ng5=[ng5 length(find(g(:) == 5))];
end

xdz=100*(0.005:.001:0.03);
figure('position',[273   138   590   658]);hold on
p1=plot(xdz,ng0,'-o','linewidth',2,'markersize',5);
p2=plot(xdz,ng1,'-o','linewidth',2,'markersize',5);
p3=plot(xdz,ng2,'-o','linewidth',2,'markersize',5);
p4=plot(xdz,ng3,'-o','linewidth',2,'markersize',5);
p5=plot(xdz,ng4,'-o','linewidth',2,'markersize',5);
p6=plot(xdz,ng5,'-o','linewidth',2,'markersize',5);
grid on;

legend([p1 p2 p3 p4 p5 p6],' g = 0',' g = 1',' g = 2',' g = 3',...
     ' g = 4',' g = 5');
 
xlabel('QC Histogram Bin Size "dz" (cm)');ylabel('Number of Lowest Goodness Values')
title(['Total Number of Survey Points: ' num2str(length(z)) ]);
