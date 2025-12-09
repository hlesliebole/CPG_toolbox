%------------------------------------------------------------
function [nsa]=UpdateSAsurvey(mpath,Mop,Survey,ndx)
%------------------------------------------------------------
%---------------------------------------------------------------------
%  Add new survey(s) to mop SA files
%---------------------------------------------------------------------

fprintf('Updating SA mop files associated with:\n%s\n',Survey(ndx).File)
%ndx=find(now - [Survey.FileDatenum] < ndays);

if ~isempty(ndx)
    
%Survey=Survey(idx); % reduce survey list to the single survey

MopStart=1; % MX border 
MopEnd=11594; % OR border

%----------------------
% now reprocess the recently modified files
%----------------------

nsurv(MopStart:MopEnd)=0;
% 

for n=ndx
    %fprintf('Survey %i of %i\n',n,size(Survey,2))
    
    [x,y,z,c,utmzone]=readSurveyFileUTM2(Survey(n).Source,Survey(n).File);
    numel(x)
  if ~isempty(x)
    ext=Survey(n).File(end-2:end); % get filename extension
    if(strcmpi(ext,'tif') == 1 || strcmpi(Survey(n).Source,'UTAir') == 1)
        % .tif and UTAir survey files already 1m averaged
    else
        % spatially average to 1m if *not* a tif file or UT air lidar
         [x,y,z,c]=SpatialAverageUTM(x,y,z,c,1);
    end
    
    [nmop]=XY2MopNumsV2(x,y,Mop);  % nearest mop area for each point
    unique(nmop)
    % save in SA struct arrays for each mop
    if ~isempty(nmop)
      for m=unique(nmop)
        if m >= MopStart && m <= MopEnd 
            nsurv(m)=nsurv(m)+1;  
            N(m).SA(nsurv(m)).Mopnum=m;
            N(m).SA(nsurv(m)).UTMzone=utmzone(1,:);
            N(m).SA(nsurv(m)).File=Survey(n).File;
            N(m).SA(nsurv(m)).FileDatenum=Survey(n).FileDatenum; 
            N(m).SA(nsurv(m)).Bytes=Survey(n).Bytes;
            N(m).SA(nsurv(m)).Source=Survey(n).Source;
            N(m).SA(nsurv(m)).Datenum=Survey(n).Datenum;
            N(m).SA(nsurv(m)).X=x(nmop == m);
            N(m).SA(nsurv(m)).Y=y(nmop == m);
            N(m).SA(nsurv(m)).Z=z(nmop == m);
            N(m).SA(nsurv(m)).Class=c(nmop == m);
        end
      end
    end
  end
end

% save SA data as individual mat files
fprintf('Saving data as SA matfiles...\n')
nsa=0; % matfile counter
for m=1:size(N,2)
  if ~isempty(N(m).SA)
    fprintf('%i\n',m)
   if exist([mpath 'M' num2str(m,'%5.5i') 'SA.mat'],'file')
    load([mpath 'M' num2str(m,'%5.5i') 'SA.mat'],'SA');
   else
    SA=[];
   end
    SA2=N(m).SA;
    SA=[SA2 SA]; % put the new data on the front end so it is the first entry for any duplicate survey
    % % sort survey list by survey date
    T=struct2table(SA);
    sortedT = sortrows(T, 'Datenum');
    SA=table2struct(sortedT)';
    
    % remove any duplicate surveys based on the filename
     [~, idx] = unique({SA.File}, 'stable');
     if numel(idx) < size(SA,2)
         fprintf(' Removing %i Duplicate Survey File Entries from SA...\n',...
             size(SA,2)-numel(idx)) 
          ddx=find(~ismember(1:size(SA,2),idx));
          for n=ddx
              fprintf('%s\n',SA(n).File)
          end
          SA=SA(idx);
     else
         fprintf('No duplicate file entries found.\n')
     end
     
    % save updated SA struct array
    save([mpath 'M' num2str(m,'%5.5i') 'SA.mat'],'SA');
    nsa=nsa+1;
  end
end
fprintf('%i matfiles updated.\n',nsa)
nsa=unique(nmop); % change nsa to unique mop numbers with data before returning

else
  fprintf('No recent Files found in SurveyMasterListWithMops.mat Survey struct array:\n%s\n',...
      fn)  
end

end
