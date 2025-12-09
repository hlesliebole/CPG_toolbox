% CheckForTrkSurveyFiles.m

%'/Volumes/group/LiDAR/VMZ2000_Truck/LiDAR_Processed_Level2/*/Beach_Only/2*ground.tif';

d=dir('/Volumes/group/topobathy/2*');

clear TrkFiles
ngps=0;

for n=1:size(d,1)
    sfile=[d(n).folder '/' d(n).name '/filtered*.*ll*.navd88'];
    s=dir(sfile);
    if(size(s,1) > 1)
        fprintf('Multiple survey files:\n')
        {s.name} 
    elseif isempty(s) % if no navd88, try adj2011 file
       sfile=[d(n).folder '/' d(n).name '/filtered*.*ll*.adj2011'];
       s=dir(sfile);
        if(size(s,1) > 1)
         fprintf('Multiple survey files:\n')
         {s.name}     
        elseif isempty(s) % if still no file, complain.
           fprintf('No survey file: %s \n',d(n).name)
        else
           %fprintf('%s\n',s.name)
           ngps=ngps+1;
           TrkFiles(ngps).folder=s.folder;
           TrkFiles(ngps).name=s.name;
        end
    else
        %fprintf('%s\n',s.name)
        ngps=ngps+1;
        TrkFiles(ngps).folder=s.folder;
        TrkFiles(ngps).name=s.name;
   end
end
