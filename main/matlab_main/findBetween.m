function [dat, inds] = findBetween(data, lw, hg)

    inds = find((data > lw) & (data < hg));
    dat = data(inds);
    
end