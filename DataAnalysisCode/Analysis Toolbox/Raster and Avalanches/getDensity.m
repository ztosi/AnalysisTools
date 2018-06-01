function [localDense, dMat, neigh] = getDensity(xys, sig, k, frs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    neigh = pdist2(xys, xys, 'euclidean', 'Smallest', k);

    xys = floor(xys);
    
    dMat = zeros(1+max(xys(:,2))-min(xys(:,2)), 1+max(xys(:,1))-min(xys(:,1)));
    
    xAx = normpdf((-2*sig):(2*sig), 0, sig);
    
    xAx = xAx' * xAx;
    [m,n] = size(dMat);
    inds = sub2ind([m,n], 1+xys(:, 2)-min(xys(:,2)), 1+xys(:,1)-min(xys(:,1)));
    
    dMat(inds) = frs./ (frs+1);
    
    dMat = conv2(dMat, xAx, 'same'); 
    
    localDense = dMat(inds);

end

