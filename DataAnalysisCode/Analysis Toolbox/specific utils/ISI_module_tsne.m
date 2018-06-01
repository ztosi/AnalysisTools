function [ mapped_states, mcolorcode, ISIs, ISI_dists, ed] = ...
    ISI_module_tsne(asdf, modules, maxISI, no_dims, pca_dims, perplexity, showfig, showhomeless)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ISIs = cellfun(@(x) nonzeros(diff(x).*(diff(x)<maxISI)), asdf(1:end-2), ...
    'UniformOutput', false);
allISIs = cell2mat(ISIs);
[~, ed] = histcounts(allISIs);
ISI_dists = cell2mat(cellfun(@(x) histcounts(x, ed, 'Normalization', ...
    'probability'), ISIs, 'UniformOutput', 0));

mcolorcode = zeros(length(asdf)-2, 1);
modszs = cellfun(@length, modules);
modules = modules(modszs>10);
modCodes = primes(100);
for ii=1:sum(modszs>10)
    mcolorcode(modules{ii}) = ii;%mcolorcode(modules{ii}) .* modCodes(ii);
end
%mcolorcode(mcolorcode==1) = 0;

if showfig
    figure;
    size(ISI_dists)
    if ~showhomeless
        mapped_states = tsne(ISI_dists(mcolorcode~=0,:), ...
            mcolorcode(mcolorcode~=0), no_dims, pca_dims, ...
            perplexity);
    else
        mapped_states = tsne(ISI_dists, ...
            mcolorcode, no_dims, pca_dims, ...
            perplexity);
    end
else
    mapped_states = tsne(ISI_dists, [], no_dims, pca_dims, ...
        perplexity);
end



end

