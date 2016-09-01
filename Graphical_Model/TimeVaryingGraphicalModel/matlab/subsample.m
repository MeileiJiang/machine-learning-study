function [ newdata ] = subsample( Data, frac )
%Generate subsample at each time point with fraction frac.
tbl = tabulate(Data.Time);
nt = ceil(tbl(:,2)*frac);
newdata = array2table(zeros(sum(nt), size(Data, 2)));
newdata.Properties.VariableNames = Data.Properties.VariableNames;
for t = 1:size(tbl, 1)
    rows = Data.Time == tbl(t,1);
    datat = Data(rows,:);
    newdata(((t-1)*nt+1):t*nt,:) = datasample(datat, nt(t));
end

end

