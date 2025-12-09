function [MopName,MopNumber]=GetMopNameNumber(MopID)

% A function to make it easy to work with both Mop character names
% based on the County, or their state-wide sequential Mop number. 
% Whichever is input to the function, it returns it and its name/number 
% counterpart as MopName and MopNumber.

County=['D    ';'OC   ';'L    ';'VE   ';'B    ';'SL   ';'MO   ';'SC   ';...
             'SM   ';'SF   ';'MA   ';'SN   ';'M    ';'HU   ';'DN   '];
CountyBeginNumber=[0 1210 1878 3063 3741 5529 6281 7202 ...
                  7530 7977 8066 8659 9172 10225 11221];

if isnumeric(MopID) % figure out MopName from the input MopNumber

    MopNumber=MopID; 
    % figure out what county the state number is in and its county number
    CountyNumbers=MopID-CountyBeginNumber;
    CountyNumbers(CountyNumbers < 1)=NaN;
    [CountyNumber,CountyIndex]=min(CountyNumbers); 
    % build MopName from its county character(s) and number
    MopName=County(CountyIndex,:);
    CountyNumeric=numel(find(MopName == ' '));
    if numel(find(MopName == ' ')) == 3
      MopName(3:5)=num2str(CountyNumber,'%3.3i');
    else
      MopName(2:5)=num2str(CountyNumber,'%4.4i');
    end

else % figure out MopNumber from the input MopName

    MopName=MopID;
    % find county that matches characters in MopName
    CountyIndex=find(strcmp(cellstr(County(:,1:5)),...
        cellstr(regexprep(MopName,'[0-9]',' '))));
    % extract the county number from MopName
    CountyNumber=str2num(regexprep(MopName,'[A-Z]',' '));
    % convert county number to state number
    MopNumber=CountyBeginNumber(CountyIndex)+CountyNumber;
    
end

end
