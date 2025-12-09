% GetVolSurveyNumber

axes(vax3);
xy=ginput(1);

 d=sqrt((Sdatenum-xy(1)).^2+(VolCutElev/1000-xy(2)).^2);
 [Y,I]=min(d);

 fprintf('Nearest Survey Date to Edit: %s\n',Sdate(I))

 CurrentSurveyNumber=find([SA(ldx).Datenum] == Sdatenum(I));