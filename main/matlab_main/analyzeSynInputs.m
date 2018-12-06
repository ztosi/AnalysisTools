function [inset, adj, adjD, mc, cs, z, zcum] = analyzeSynInputs(n, synConMat)


inset = (synConMat(:,n) ~=0) & (synConMat(n,:)' == 0);
adjD = double(synConMat(inset, inset) ~= 0);
adj = (adjD' + adjD)~=0;
size(adj)
nnz(adj)/(size(adj,1)*(size(adj,2)-1))
mc = maximalCliques(adj~=0);
[cs, z, zcum] = findFlows(mc, adjD);

end