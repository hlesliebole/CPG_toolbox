% find and plot all the July 2016 and Jul 2024 data
clearvars 
close all

load ProfileWidthsVolumes.mat
MopNum=530:780;
figure('position',[ 144          62        1192         735]);

for cy=2023:2024 % comparison year
    nm=0;
for m=530:780

nm=nm+1;

bw1(nm)=NaN;
bv1(nm)=NaN;
idx=find(month(Zt(nm).mon) == 7 & year(Zt(nm).mon) == cy);
if isempty(idx)
idx=find(month(Zt(nm).mon) == 6 & year(Zt(nm).mon) == cy);
end
if isempty(idx)
idx=find(month(Zt(nm).mon) == 8 & year(Zt(nm).mon) == cy);
end
 
if ~isempty(idx)
 bw1(nm)=BWmax(nm).mon(idx);
 bv1(nm)=BVmax(nm).mon(idx);
end

% 2024

bw2(nm)=NaN;
bv2(nm)=NaN;
idx=find(month(Zt(nm).mon) == 7 & year(Zt(nm).mon) == cy+1);
if isempty(idx)
idx=find(month(Zt(nm).mon) == 6 & year(Zt(nm).mon) == cy+1);
end
if isempty(idx)
idx=find(month(Zt(nm).mon) == 8 & year(Zt(nm).mon) == cy+1);
end
 
if ~isempty(idx)
 bw2(nm)=BWmax(nm).mon(idx);
 bv2(nm)=BVmax(nm).mon(idx);
end


end

%figure;
subplot(4,1,1);grid on;hold on;
plot(MopNum,bw2,'.-');title('Jul 2020 vs Jul 2024 beach width');
set(gca,'xlim',[MopNum(1) MopNum(end)]);grid on;
subplot(4,1,2);grid on;hold on;
plot(MopNum,bw2-bw1,'.-');title('Jul 2024 vs Jul 2020 beach width change');
set(gca,'xlim',[MopNum(1) MopNum(end)]);grid on;

%figure;
subplot(4,1,3);
plot(MopNum,bv2,'.-');title('Jul 2020 vs Jul 2024 beach width');
set(gca,'xlim',[MopNum(1) MopNum(end)]);grid on;hold on;
subplot(4,1,4);grid on;hold on;
plot(MopNum,bv2-bv1,'.-');title('Jul 2020 vs Jul 2024 beach width change');
set(gca,'xlim',[MopNum(1) MopNum(end)]);grid on;

end

