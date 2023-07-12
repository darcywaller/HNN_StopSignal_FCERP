function [p,actualdiff] = mc_ttest(sample1,sample2,mciter)

actualdiff = mean(sample1) - mean(sample2);
datavector = [sample1; sample2];
n1 = length(sample1);
n2 = length(sample2);
for mci = 1:mciter-1
    reorder = datavector(randperm(length(datavector)));
    sample1 = reorder(1:n1);
    sample2 = reorder(n1+1:end);
    mcdiff(mci) = mean(sample1) - mean(sample2);
end
mcdiff = sort([mcdiff actualdiff]);
if actualdiff > 0
    p = 2*(1-find(mcdiff == actualdiff,1,'first')/mciter);
else p = 2*(find(mcdiff == actualdiff,1,'first')/mciter);
end
if p > 1; p = 1; end