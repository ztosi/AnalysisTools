function [totFlow, threshFlow, cumuFlow] = findFlows(cliqueSet, adjD)


totFlow = zeros(size(cliqueSet,2),1);
threshFlow = zeros(size(cliqueSet))';
cumuFlow = zeros(size(threshFlow));
for ii=1:size(cliqueSet,2)
    indx = find(cliqueSet(:,ii));
    threshFlow(ii,:) = cliqueSet(:,ii);
    for jj=1:100
        threshFlow(ii,indx) = threshFlow(ii,indx) * adjD(indx, indx);
        cumuFlow(ii,indx) = cumuFlow(ii,indx) + threshFlow(ii,indx);
        threshFlow(ii, threshFlow(ii,:)>1) = 1;
        totFlow(ii) = totFlow(ii) + sum(threshFlow(ii,indx));
    end
end