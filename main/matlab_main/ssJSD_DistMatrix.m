function [ jsdMat ] = ssJSD_DistMatrix ( asdf, kernel )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = length(asdf)-2;

intspkints = cellfun(@(x) diff(x)', ...
    asdf(1:end-2), 'UniformOutput', 0);

jsdMat = zeros(N);

for ii=1:N
    disp(ii)
    parfor jj = (ii+1):N
        jsdMat(ii,jj) = findMaxDistSimilarity(intspkints{ii}, ...
            intspkints{jj}, 'Kernel', kernel, 'Reshape', 'Cutoff', 50000);
%          jsdMat(ii,jj) = findMaxDistSimilarity(intspkints{ii}, ...
%              intspkints{jj}, 'Kernel', kernel, 'Cutoff', 50000);
    end
end

jsdMat = jsdMat + jsdMat';

end

