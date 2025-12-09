% first find all the official jumbo survey data indexes
jdx=find(contains({SA.File},'umbo'));

% loop through those and if there are any Trk or AtvMR surveys within 2 days
for Surv=jdx 
 ldx=find((strcmp({SA.Source},'AtvMR') | strcmp({SA.Source},'Trk')) &...
 ([SA.Datenum] >= SA(Surv).Datenum-2 & [SA.Datenum] <= SA(Surv).Datenum+2));
%  add any LiDAR x,y,z data to the jumbo data fields
 if ~isempty(ldx)
  for SameSurv=ldx
      SA(Surv).X=vertcat(SA(SameSurv).X,SA(Surv).X);
      SA(Surv).Y=vertcat(SA(SameSurv).Y,SA(Surv).Y);
      SA(Surv).Z=vertcat(SA(SameSurv).Z,SA(Surv).Z);
  end
 end
end

% now reduce to just (LiDAR enhanced) jumbos
SAR=SA(jdx);
