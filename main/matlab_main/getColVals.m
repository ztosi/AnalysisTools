function vec = getColVals(dat, colV)

    vec = zeros(size(dat,1), 1);
    vals = unique(colV(:));
    for ii=1:length(vals)
       ind = vals(ii);
       wh = find(colV == ind);
       vec(wh) = dat(wh, ind);
    end

end