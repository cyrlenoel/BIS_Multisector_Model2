function [results] = cell2List(cc)
for c =1:length(cc)
    if c == 1
        results = cc{c};
    else
        results = strcat(results,',',cc{c});
    end
end