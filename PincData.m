load M00589SA.mat
n=size(SA,2);
for m=n-10:n
    fprintf('%s  %s %s\n',datestr(SA(m).Datenum),SA(m).Source,SA(m).File)
end
    