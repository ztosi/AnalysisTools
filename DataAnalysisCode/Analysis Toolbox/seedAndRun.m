function [acts, seeds] = seedAndRun(wtMat, posSeeds, nseed, iters, runs, thresh)


acts = zeros(size(wtMat,1), iters+1);
seeds = zeros(nseed, runs);


for ii=1:runs
    vec = zeros(1, size(wtMat,1));
    timer = zeros(1, size(wtMat,1));
    seeds(:, ii) = posSeeds(randperm(length(posSeeds), nseed));
    vec(seeds(:, ii)) = 1;
    acts(:, 1) = acts(:,1) + vec';
    timer(vec == 1) = iters;
    for jj=2:iters+1
        vec = vec * wtMat;
        vec = (vec > thresh) & (timer<=0);
        acts(:, jj) = acts(:,jj) + vec';
        timer = timer - 1;
        timer(vec==1) = iters;
    end
    
end

acts = acts ./ runs;

end

