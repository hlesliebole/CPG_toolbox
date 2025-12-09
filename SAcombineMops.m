function CS=SAcombineMops(MopStart,MopEnd)

% combines gridded SG info in a reach of Mops and returns the combined grid
% data in a combined CG struct array of utm X,Y and Z variables.  

% initialize combined struct array with the SG struct array for the first 
% Mop.

CS=[];

% load remaining mop SG arrays and add data to existing file 
%  entries in CG or add a new entry to it.
for MopNumber=MopStart:MopEnd
     
 matfile=['M' num2str(MopNumber,'%5.5i') 'SA.mat' ];

 if exist(matfile,'file')
  load(matfile,'SA');
  fprintf('Loaded %s with %i surveys\n',matfile,size(SA,2));
 
  if isempty(CS)
      
      CS=SA;
  else
      
    for n=1:size(SA,2)
      SameSurv=find(strcmp({CS.File},SA(n).File));
      if ~isempty(SameSurv)
          CS(SameSurv).Mopnum=vertcat(CS(SameSurv).Mopnum,SA(n).Mopnum);
          CS(SameSurv).X=vertcat(CS(SameSurv).X,SA(n).X);
          CS(SameSurv).Y=vertcat(CS(SameSurv).Y,SA(n).Y);
          CS(SameSurv).Z=vertcat(CS(SameSurv).Z,SA(n).Z);
          CS(SameSurv).Class=vertcat(CS(SameSurv).Class,SA(n).Class);
      else
          CS(size(CS,2)+1)=SA(n);
      end
    end
    
  end
  
 else
   fprintf('%s does not exist.\n',matfile);   
 end
end

% sort CS by survey date
    if size(CS,2) > 1
     T=struct2table(CS);
     sortedT = sortrows(T, 'Datenum');
     CS=table2struct(sortedT)';
    end
    
end

