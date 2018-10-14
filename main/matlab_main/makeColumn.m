function [out] = makeColumn(input)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if (iscell(input))
    out = cellfun(@(x) x(:), input, 'UniformOutput', 0);
else
    out = input(:);
end

end

