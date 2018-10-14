function [ totSsts, inSumsts, outSumsts, allsynssts ] = syncumsum( wtMat )


if iscell(wtMat)
    numSets = size(wtMat,1);
    sets = cell(numSets,4);
    for i=1:numSets
        [sets{i,1}, sets{i,2}, sets{i,3}, sets{i,4}] = ...
            cumsum50(wtMat{i});
    end    
else
    numSets = size(wtMat,3);
    sets = cell(numSets,4);
    for i=1:numSets
        [sets{i,1}, sets{i,2}, sets{i,3}, sets{i,4}] = ...
            cumsum50(wtMat(:,:,i));
    end
end

if numSets == 1
    base = sets(1, 1:4);
    sets = cell(1000, 4);
    if iscell(wtMat)
        for i=1:1000
            ss = randperm(size(wtMat{1},1), uint32(size(wtMat{1},1)/10));
            [sets{i,1}, sets{i,2}, sets{i,3}, sets{i,4}] = ...
                cumsum50(wtMat{1}(ss,ss));
        end
    else
        for i=1:1000
            ss = randperm(size(wtMat,1), uint32(size(wtMat,1)/10));
            [sets{i,1}, sets{i,2}, sets{i,3}, sets{i,4}] = ...
                cumsum50(wtMat(ss,ss));
        end    
    end
    nn = length(base{1});
    ns = length(base{4});
    totSsts = zeros(3, nn);
    inSumsts = zeros(3, nn);
    outSumsts = zeros(3, nn);
    allsynssts = zeros(3, ns);
    
    totSsts(1,:) = base{1,1};
    inSumsts(1,:) = base{1,2};
    outSumsts(1,:) = base{1,3};
    allsynssts(1,:) = base{1,4};
    
    
    totl = zeros(1000, nn);
    inl = zeros(1000, nn);
    oul = zeros(1000, nn);
    
    for i=1:1000
       totl(i,:) = stretchWithCopies(sets{i,1}, nn);
       inl(i,:) = stretchWithCopies(sets{i,2}, nn);
       oul(i,:) = stretchWithCopies(sets{i,3}, nn);
    end
    [totl, totSsts(2,:)] = semistd(totl);
    [inl, inSumsts(2,:)] = semistd(inl);
    [oul, outSumsts(2,:)] = semistd(oul);
    
    totSsts(3,:) = totl;
    inSumsts(3,:) = inl;
    outSumsts(3,:) = oul;
    
    mxa = 10*max(cellfun(@numel, sets(:,4)));
    mxc = cell(1000,1);
    [mxc{:,1}] = deal(mxa);
    alls = cellfun(@stretchWithCopies, sets(:,4), mxc, 'UniformOutput', ...
        false);
    
    [alls, allu] = semistd(cell2mat(alls));
    
    allsynssts(2,:) = stretchWithCopies(alls, ns);
    allsynssts(3,:) = stretchWithCopies(allu, ns);  
else
    
    
    
    
end




end

