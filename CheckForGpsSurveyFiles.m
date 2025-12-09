% CheckForGpsSurveyFiles.m

d=dir('/Volumes/group/topobathy/2*');

clear GpsFiles
ngps=0;

for n=1:size(d,1)
    sfile=fullfile(d(n).folder,d(n).name,'filtered*.*ll*.navd88');
    s=dir(sfile);
    if(size(s,1) > 1)
        fprintf('Multiple survey files:\n')
        {s.name} 
    elseif isempty(s) % if no navd88, try adj2011 file
       sfile=fullfile(d(n).folder,d(n).name,'filtered*.*ll*.adj2011');
       s=dir(sfile);
        if(size(s,1) > 1)
         fprintf('Multiple survey files:\n')
         {s.name}     
        elseif isempty(s) % if still no file, complain.
           fprintf('No survey file: %s \n',d(n).name)
        else
           %fprintf('%s\n',s.name)
           ngps=ngps+1;
           GpsFiles(ngps).folder=s.folder;
           GpsFiles(ngps).name=s.name;
        end
    else
        %fprintf('%s\n',s.name)
        ngps=ngps+1;
        GpsFiles(ngps).folder=s.folder;
        GpsFiles(ngps).name=s.name;
   end
end
