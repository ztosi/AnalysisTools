function [mnSyn, totSyn, mnUnq, totUnq, mnRed, totRed, n] = findSynVals(r_n, t_n, pilData)

    recInds = find(pilData(:,1) == r_n);
    [transInd, transNo] = find(pilData(:,2:3) == t_n);

    [inds, ~, tind] = intersect(recInds, transInd);
    transNo = transNo(tind);
    
    pilData = pilData(inds,:);
    for ii=1:length(inds)
        if transNo(ii) == 3
           tmpT = pilData(ii,3);
           pilData(ii,3) = pilData(ii,2);
           pilData(ii,2) = tmpT;
           tmpUI = pilData(ii,7);
           pilData(ii,7) = pilData(ii,6);
           pilData(ii,6) = tmpUI;
        end
    end
    
    mnSyn = mean(pilData(:,4));
    totSyn = sum(pilData(:,4));
    totUnq = sum(pilData(:,6));
    mnUnq = mean(pilData(:,6));
    totRed = sum(pilData(:,5));
    mnRed = mean(pilData(:,5));
    n = length(inds);


end