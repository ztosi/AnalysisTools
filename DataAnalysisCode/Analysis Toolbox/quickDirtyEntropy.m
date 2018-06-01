function ent = quickDirtyEntropy(data)

    [bh, ~] = histcounts(data, 'Normalization', 'Probability');
    nz=bh~=0;
    ent = -sum(bh(nz) .* log2(bh(nz)));
    
end