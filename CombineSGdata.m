function CG=CombineSGdata(mpath,MopStart,MopEnd)

% combines gridded SG info in a reach of Mops and returns the combined grid
% data in a combined CG struct array of the utm X,Y and Z variables. 

% mpath = full path name to the cgpmop M*SG.mat file folder
%   eg. on a mac /volumes/group/MOPS/

% The CG struct array has the same fields as the SG array, but the
%  CG(N).Mopnum field now lists all the combined mop numbers that 
%  have data for survey N rather than a single mop number.

% initialize combined struct array with the SG struct array for the first 
% Mop.

CG=[];

% load remaining mop SG arrays and add data to existing file 
%  entries in CG or add a new entry to it.
for MopNumber=MopStart:MopEnd
     
 matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SG.mat' ];

 if exist(matfile,'file')
  load(matfile,'SG');
  fprintf('Loaded %s with %i surveys\n',matfile,size(SG,2));
 
  if isempty(CG)
      
      CG=SG;
  else
      
    for n=1:size(SG,2)
      SameSurv=find(strcmp({CG.File},SG(n).File));
      if ~isempty(SameSurv)
          CG(SameSurv).Mopnum=vertcat(CG(SameSurv).Mopnum,SG(n).Mopnum);
          CG(SameSurv).X=vertcat(CG(SameSurv).X,SG(n).X);
          CG(SameSurv).Y=vertcat(CG(SameSurv).Y,SG(n).Y);
          CG(SameSurv).Z=vertcat(CG(SameSurv).Z,SG(n).Z);
          CG(SameSurv).Class=vertcat(CG(SameSurv).Class,SG(n).Class);
      else
          CG(size(CG,2)+1)=SG(n);
      end
    end
    
  end
  
 else
   fprintf('%s does not exist.\n',matfile);   
 end
end

% sort CG by survey date
    if size(CG,2) > 1
     T=struct2table(CG);
     sortedT = sortrows(T, 'Datenum');
     CG=table2struct(sortedT)';
    end
    
end

