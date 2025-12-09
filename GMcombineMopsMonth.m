function CG=GMcombineMopsMonth(MopStart,MopEnd,Mon)

% combines gridded SG info in a reach of Mops and returns the combined grid
% data in a combined CG struct array of the utm X,Y and Z variables. 

% The CG struct array has the same fields as the SG array, but the
%  CG(N).Mopnum field now lists all the combined mop numbers that 
%  have data for survey N rather than a single mop number.

% initialize combined struct array with the SG struct array for the first 
% Mop.

CG=[];

% load remaining mop SG arrays and add data to existing file 
%  entries in CG or add a new entry to it.
for MopNumber=MopStart:MopEnd
     
 matfile=['M' num2str(MopNumber,'%5.5i') 'GM.mat' ];

 if exist(matfile,'file')
  load(matfile,'GM');
  fprintf('Loaded %s \n',matfile);
  
 
  if isempty(CG)
      CG.X2D=GM.MM(Mon).X2D;
      CG.Y2D=GM.MM(Mon).Y2D;
      CG.Z2Dmean=GM.MM(Mon).Z2Dmean;
      CG.Z2Dmedian=GM.MM(Mon).Z2Dmedian;
  else
      CG.X2D=vertcat(CG.X2D,GM.MM(Mon).X2D);
      CG.Y2D=vertcat(CG.Y2D,GM.MM(Mon).Y2D);
      CG.Z2Dmean=vertcat(CG.Z2Dmean,GM.MM(Mon).Z2Dmean);
      CG.Z2Dmedian=vertcat(CG.Z2Dmedian,GM.MM(Mon).Z2Dmedian);
  end

 end
end

end

