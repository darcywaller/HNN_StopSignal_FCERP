% Grubbs' Test for Outliers
% by Jan R. Wessel 2009
% Wessel@nf.mpg.de
%
function [ val , idents , critical ] = outlier_grubbs( input , alpha , sided )

if nargin == 2
    sided = 'both';
elseif nargin == 1
    sided = 'both';
    alpha = 0.05;
end

n_grubbs = length(input);

% Critical Value
if strcmpi(sided,'both')
    multiplier = 2;
else multiplier = 1;
end
grubbs_t = (tinv(alpha/(multiplier*n_grubbs),(n_grubbs - 2)))^2;
g2denom = (n_grubbs - 2) + grubbs_t;
g2 = sqrt(grubbs_t / g2denom);
g1 = (n_grubbs - 1) / sqrt(n_grubbs);
grubbs = g2*g1;

% Test statistic
outcnt = 0; val = []; id = []; originput = input; idents = [];
if strcmpi(sided,'both')   % Two-sided
    gvalues = abs( (input - mean(input))/std(input) );
    while max(gvalues) > grubbs
        outcnt = outcnt + 1;
        if length(find(gvalues == max(gvalues))) > 1
            temp = (find(gvalues == max(gvalues)));
            id(outcnt) = temp(1);
        else id(outcnt) = (find(gvalues == max(gvalues)));
        end
        val(outcnt) = input(id(outcnt));
        input = setxor(originput,val);
        gvalues = ( (input - mean(input))/std(input) );
    end
elseif strcmpi(sided,'positive')  % One-Sided, positive
    gvalues = ( (input - mean(input))/std(input) );
    while max(gvalues) > grubbs
        outcnt = outcnt + 1;
        if length(find(gvalues == max(gvalues))) > 1
            temp = (find(gvalues == max(gvalues)));
            id(outcnt) = temp(1);
        else id(outcnt) = (find(gvalues == max(gvalues)));
        end
        val(outcnt) = input(id(outcnt));
        input = setxor(originput,val);
        gvalues = ( (input - mean(input))/std(input) );
    end
elseif strcmpi(sided,'negative')   % One-Sided, negative
    gvalues = ( (input - mean(input))/std(input) );
    while min(gvalues) < -grubbs
        outcnt = outcnt + 1;
        if length(find(gvalues == max(gvalues))) > 1
            temp = (find(gvalues == max(gvalues)));
            id(outcnt) = temp(1);
        else id(outcnt) = (find(gvalues == max(gvalues)));
        end
        val(outcnt) = input(id(outcnt));
        input = setxor(originput,val);
        gvalues = ( (input - mean(input))/std(input) );
    end
else
    error('Direction method not recognized! Use "positive", "negative", or "both" as stringinputs to outlier_grubbs.')
end

for i = 1:numel(val)
    idents(i) = find(originput == val(i));
end

critical = grubbs;
if isempty(val)
    disp('The sample does not contain outliers.');
end