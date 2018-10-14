function [actual, nulls] = inputConnectivity(wtMat, noNulls)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
adj=wtMat~=0;
kIn = sum(adj);

[N,~] = size(wtMat);
nulls = zeros(N, noNulls);
actual = zeros(N,1);
parfor jj=1:noNulls
    for ii=1:N
        samp = randperm(N, kIn(ii));
        nulls(ii,jj) = nnz(adj(samp, samp)) / (kIn(ii)*(kIn(ii)-1));
    end
end

for ii=1:N
    samp = adj(:,ii);
   actual(ii) = nnz(adj(samp, samp)) / (kIn(ii)*(kIn(ii)-1));
end

end

