function [ent_tot, ent_pc, ent_d] = clusterEntropy(x, idx, centers, kappa)

    ent_pc = zeros(length(centers),1);
    N = size(x, 1);
    K = size(centers,1);
    D = pdist2(x,x);
    D = exp(-(D./kappa).^2);
    D = bsxfun(@rdivide, D, sum(D));
    ent_d = sum(-sum(D .* log2(D)))/N;
    ent_tot=0;
    for ii=1:K
       cdists = bsxfun(@minus, x, centers(ii,:));
       cdists = sqrt(sum(cdists.^2,2)); 
       cdists = exp(-(cdists./kappa).^2);
       cdists = cdists./sum(cdists);
       ent_pc(ii) = -sum(cdists .* log2(cdists)) ;%+ ... 
           %(sum(idx==ii)/N .* log2(sum(idx==ii)/N) );
       crs = bsxfun(@rdivide, D(:, idx==ii), cdists);
       ent = sum(sum(D(:, idx==ii) .* log2(crs)));
       
       ent_tot = ent_tot + (ent + ent_pc(ii))/N  + ... 
           (sum(idx==ii)/N .* log2(sum(idx==ii)/N));
    end


end