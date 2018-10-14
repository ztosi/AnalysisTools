function [cent, ent] = cumEntropy(mat)
[N, M] = size(mat);
mat = sort(mat);
mat2 = bsxfun(@plus, mat, -min(mat));
mx = max(mat2);
mat2 = bsxfun(@rdivide, mat2, mx);
dfs = diff(mat2,1);
% mat = cumsum(mat);
% mx = max(mat);
% mat = bsxfun(@rdivide, mat, mx);

mat2 = 1/N:1/N:(1-(1/N));
mat2 = 1-repmat(mat2', 1, M); 
figure; plot(1-cumsum(dfs(:,1)), mat2(:,1));
figure; plot(mat(2:end,1)', (1/N)./diff(mat(:,1)));
sum((1/N)./diff(mat(:,1)))
cent = -sum(mat2 .* dfs .*log2(mat2) );

end