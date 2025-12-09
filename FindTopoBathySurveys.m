% Build SM  Morpho struct files from SG grid files

%mpath='/volumes/group/MOPS/';
DefineMopPath

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');

% input Mop number

% Mop1=571;
% Mop2=594;
% Mop1=578;
% Mop2=584;
% Mop1=550;
% Mop2=555;
Mop1=1;%2;%550;%229; % MX border 
Mop2=11594; % OR border
%Mop2=500;
% Mop1=578;Mop2=584;% BN+ruby2d
%Mop1=664;Mop2=683;

 for MopNumber=Mop1:Mop2
  SM=[];
  ntp=0;  
    
    
    matfile=[mpath 'M' num2str(MopNumber,'%5.5i') 'SM.mat' ];

    if exist(matfile,'file')
     load(matfile,'SM');
    if size(SM,2) > 0
   
      
    for n=1:size(SM,2)
        if SM(n).Datenum >= datenum(2016,1,1)
        if size(SM(n).Z1Dmean,1) > 0
        if min(SM(n).Z1Dmean) < -6 && max(SM(n).Z1Dmean) > 1
            if contains(SM(n).File,'umbo') || strcmp(SM(n).Source,'USACE')...
                || strcmp(SM(n).Source,'CCC' )
              ntp=ntp+1;
            end
        end
        end
        end
    end
    
    
    end
    end
    NumTopoBathy(MopNumber)=ntp;
    fprintf('Mop %i has %i topobathy surveys\n',MopNumber,ntp)
 end
  
