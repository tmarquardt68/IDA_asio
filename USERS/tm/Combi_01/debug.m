n=0;
for q2 = 1:16
    for q = 1:5
        n=n+1;
        eval(['t(' num2str(n) ')=raw_' num2str(q2) '_' num2str(q) '.timeStamp;']);
    end
end
