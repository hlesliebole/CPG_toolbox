% makes a SHOALS survey roughness struct array for a range of Mop areas.

RC=[];
%m=755:759;% grandview = 757ish
%m=474:495;% LJ Childrens Pool to South End LJ Shores Beach 
%m=660:723;% Cardiff SB to Moonlight SB
%m=637:667;%solana beach
m=570:600;% North Torrey
for MopNumber=m
 
 fprintf('Mop %i\n',MopNumber)    
 RS=GetMopShoalsRoughness(MopNumber);
 fprintf('RS has %i Roughness points.\n',numel([RS.Sigma]))  
 if isempty(RC)
     RC=RS;
 else
     RC(MopNumber-min(m)+1)=RS;
 end
 
end

matfile=['M' num2str(m(1),'%5.5i') 'to' num2str(m(end),'%5.5i') 'RC.mat'];
fprintf('Saving RC struct array to %s\n',matfile)  
save(matfile,'RC');
