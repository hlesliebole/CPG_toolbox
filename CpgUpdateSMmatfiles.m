
clearvars
CpgDefineMopPath

load SurveyMasterListWithMops.mat % Load current Survey info
load('MopTableUTM.mat','Mop');  % Load "Mop" table array

% find all the Mop numbers with survey data
MopNums=unique([Survey.NearestMops]);
MopNums=sort(MopNums);

for MopNumber=MopNums
    m=MopNumber;
    % find surveys with data for this MopID
ndx=find(cellfun(@(v)any(v(:) == m),{Survey.NearestMops}));
 
% read in the SM file if it exists

SMfiles=[];
if exist([ mpath 'M'  num2str( m , '%5.5i' )  'SM.mat' ],'file')  
  load([ mpath 'M'  num2str( m , '%5.5i' )  'SM.mat' ],'SMfiles');
else
  SMfiles=[];
  SM=[];
end

   if size(SMfiles,2) > 0
   % combine filedatenums and filenames into single strings for comparison
   % of both the file names and creation dates at the same time
   Mdf=strcat(string([SMfiles.FileDatenum]),{SMfiles.File}); % SG string
   Sdf=strcat(string([Survey(ndx).FileDatenum]),{Survey(ndx).File}); % Survey string
   
   ndel=find(~ismember(Mdf,Sdf)); % SM file that is not a member of Survey
   nadd=find(~ismember(Sdf,Mdf)); % Survey file not a member of SM
   if numel(ndel) > 0
       fprintf(' Mop %i Deleting %i surveys\n',m,numel(ndel))
   end
   if numel(nadd) > 0
       fprintf(' Mop %i Adding %i surveys\n',m,numel(nadd))
   end
   
   else
       fprintf('No SM data for Mop %i.  %i surveys being added.\n',m,numel(ndx))
       ndel=[];
       nadd=ndx;
   end

% load existing SM struct array 
if numel(ndel) > 0 || numel(nadd) > 0
 if exist([ mpath 'M'  num2str( m , '%5.5i' )  'SM.mat' ],'file')  
  load([ mpath 'M'  num2str( m , '%5.5i' )  'SM.mat' ],'SM');
 else
  SM=[];
 end
end

% delete changed surveys
if numel(ndel) > 0
    SM(ndel)=[];
    SMfiles=struct('File',{SM.File},'FileDatenum',...
            num2cell([SM.FileDatenum]));
    save([ mpath 'M'  num2str( m , '%5.5i' )  'SM.mat' ],'SM','SMfiles');
end
   
if numel(nadd) > 0
    CpgUpdateSM(MopNumber);
end

end
