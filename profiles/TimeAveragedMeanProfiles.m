function [TZdatetime,TZ1Dmon,TZ1Dqtr,TZ1Dsea,TZ1Dann,TZ1Dglo]=...
                     TimeAveragedMeanProfiles(Zdatetime,Z1D)

% Calculates progressive monthly > quarterly > 2 season > annual > global 
%   time averaged profiles from the input Z1D profile matrix

% input
%
% Zdatetime(M) = matlab datetimes of the M surveys in Z1D
% Z1D(M,N) = M survey profiles, with N cross-shore points

%
% returns 
%
%  TZdatetime = struct array of timescale-centered datetimes of averaged
%                profiles. eg. TZdatetime.mon has datetimes for the middle
%                of the months with time avergaed data.
%
%  TZ1Dmon(A,N) = monthly mean profiles for the A months with data 
%  TZ1Dqtr(B,N) = quarterly (OND,JFM,AMJ,JAS) mean profiles for the 
%                   B quarters with data (derived from monthly means)
%  TZ1Dsea(C,N) = 2 season (ONDJFM,AMJJAS) mean profiles for the C seasons 
%                   with data (derived from quarterly means)
%  TZ1Dann(D,N) = annual (Oct [yr-1] to Sep yr) mean profiles for the D years 
%                 with data (derived from season means, NaN if one season missing)  
%  TZ1Dglo(1,N) = global mean profile (derived from annual means)

% survey datetimes
Zdt=Zdatetime;

%% year-month means
n=0;
for y=year(Zdt(1)):year(Zdt(end))
   for m=1:12
       idx=find(year(Zdt) == y & month(Zdt) == m);
       if numel(idx) > 0
           n=n+1;
           Zymdate(n)=datetime(y,m,15,0,0,0);
           
           if numel(idx) == 1
               TZ1Dmon(n,:)=Z1D(idx,:);
           elseif numel(idx) > 1
               TZ1Dmon(n,:)=mean(Z1D(idx,:),'omitnan');
           end
       end
   end
end

TZdatetime.mon=Zymdate;

%% year-quarter means
Zdt=TZdatetime.mon;
n=0;
for y=year(Zdt(1)):year(Zdt(end))
   for q=1:4
       idx=find(year(Zdt) == y & ceil(month(Zdt)/3) == q);
       if numel(idx) > 0
           n=n+1;
           Zyqdate(n)=datetime(y,3*(q-1)+2,15,0,0,0);
           
           if numel(idx) == 1
               TZ1Dqtr(n,:)=TZ1Dmon(idx,:);
           elseif numel(idx) > 1
               TZ1Dqtr(n,:)=mean(TZ1Dmon(idx,:),'omitnan');
           end
       end
   end
end

TZdatetime.qtr=Zyqdate;

%% year- 2 season means (oct-sep beach years)
Zdt=TZdatetime.qtr+calmonths(3); % advance time by 3 months to shift logic to beach years
n=0;
for y=year(Zdt(1)):year(Zdt(end))
   for s=1:2
       idx=find(year(Zdt) == y & ceil(month(Zdt)/6) == s);
       if numel(idx) > 0
           n=n+1;
           % mid jan 1 and mid jul 1 are 2 season centers
           Zysdate(n)=datetime(y,6*(s-1)+1,1,0,0,0);
           
           if numel(idx) == 1
               TZ1Dsea(n,:)=TZ1Dqtr(idx,:);
           elseif numel(idx) > 1
               TZ1Dsea(n,:)=mean(TZ1Dqtr(idx,:),'omitnan');
           end
       end
   end
end

TZdatetime.sea=Zysdate;

%% annual means (oct-sep beach years)
Zdt=TZdatetime.sea;

Zydate=datetime([],[],[]);
TZ1Dann=[];
n=0;
for y=year(Zdt(1)):year(Zdt(end))
       idx=find(year(Zdt) == y);
       % need both seasons to make an annual value
       if numel(idx) == 2
           n=n+1;
           % beach year time center is Apr 1
           Zydate(n)=datetime(y,4,1,0,0,0);
           % dont omitnan for annual avg. want both seasons to have a
           %   vaild xshore point
           TZ1Dann(n,:)=mean(TZ1Dsea(idx,:));
       end
end

if ~isempty(Zydate)
TZdatetime.ann=Zydate;
end
%% global mean

if exist('TZ1Dann','var')
TZ1Dglo=mean(TZ1Dann,'omitnan');
end

% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% % monthly means across all years
% n=0;
% for m=1:12
%        n=n+1;
%        Zm(n,:)=Zym(1,:)*NaN;
%        idx=find(Zmonth == m);
%        if numel(idx) == 1
%            Zm(n,:)=Zym(idx,:);
%        elseif numel(idx) > 1
%            Zm(n,:)=mean(Zym(idx,:),'omitnan');
%        end
% end
% 
% 
% % 3-month running means across all years
% n=0;
% for m3=1:12
%        n=n+1;
%        Zm3(n,:)=Zym(1,:)*NaN;
%        if m3 == 1
%            m=[1 2 12];
%        elseif m3 == 12
%            m=[1 11 12];
%        else
%            m=m3-1:m3+1;
%        end
% 
%        idx=find(Zmonth == m(1) | Zmonth == m(2) | Zmonth == m(3));
%        if numel(idx) == 1
%            Zm3(n,:)=Zym(idx,:);
%        elseif numel(idx) > 1
%            Zm3(n,:)=mean(Zym(idx,:),'omitnan');
%        end
% end
% 
% 
% % quarterly means from monthly means
% n=0;
% for q=1:4
%        n=n+1;
%        Zq(n,:)=Zm(1,:)*NaN;      
%        idx=(q-1)*3+1:q*3;
%        Zq(n,:)=mean(Zm(idx,:),'omitnan');
% end
% 
% % figure;
% % hold on;
% % for m=1:4
% %     plot(X1Dt,Zq(m,:),'-','linewidth',2)
% %     hold on;
% % end
% % legend
% 
% % seasonal means from quarterly means
% Zs(1,:)=mean(Zq(1:2,:),'omitnan');
% Zs(2,:)=mean(Zq(3:4,:),'omitnan');
% 
% % figure;
% % hold on;
% % for m=1:2
% %     plot(X1Dt,Zs(m,:),'-','linewidth',2)
% %     hold on;
% % end
% % legend
% 
% % global mean
% Zg=mean(Zs,'omitnan');


end