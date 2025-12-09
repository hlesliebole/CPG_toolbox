SGf=dir('M*SG.mat');

for n=1:size(SGf,1)

matfile=fullfile(SGf(n).folder,SGf(n).name)
load(matfile,'SG')

% 
for s=1:size(SG,2)
  SG(s).Class=zeros(size(SG(s).X));
end

  eval(['save ' matfile ' SG']);
end
