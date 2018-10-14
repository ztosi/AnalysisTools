function [stats] = interModStats( mods, wtMat, ei )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    stats1 = zeros(length(mods));
    for ii=1:length(mods)
        for jj = 1:length(mods)
            stats1(ii,jj) = sum(sum(abs(wtMat(intersect(mods{ii}, find(ei)),...
                mods{jj}))));
        end
    end
    stats = {stats1};
    stats{2} = bsxfun(@rdivide, stats1, sum(stats1));
    stats{3} = bsxfun(@rdivide, stats1, cellfun(@length, mods));


end

