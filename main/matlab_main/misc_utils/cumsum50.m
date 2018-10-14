function [ totSum, inSum, outSum, allsyns] = cumsum50( wtMat )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    wtMat(isnan(wtMat)) = 0;
    wtMat = abs(wtMat);
    inSum = sum(wtMat);
    outSum = sum(wtMat, 2)';
    totSum = inSum + outSum;
    
    allsyns = sort(nonzeros(wtMat(:)), 'descend');
    [inSum, inI] = sort(inSum, 'descend');
    [outSum, outI] = sort(outSum, 'descend');
    [totSum, totI] = sort(totSum, 'descend');

    su = sum(inSum); % doesn't matter...
    
    inSum = cumsum(inSum)./sum(inSum);
    outSum = cumsum(outSum)./sum(outSum);
    totSum = cumsum(totSum)./(sum(totSum));
    allsyns = cumsum(allsyns)./sum(allsyns);

%     [~,t50] = min(abs(totSum-.5));
%     [~,i50] = min(abs(inSum-.5));
%     [~,o50] = min(abs(outSum-.5));
%     [~,a50] = min(abs(allsyns-.5));
    
%     [~,t70] = min(abs(totSum-.7));
%     [~,i70] = min(abs(inSum-.7));
%     [~,o70] = min(abs(outSum-.7));
%     [~,a70] = min(abs(allsyns-.7));
% 
%     N = size(wtMat,1);
% 
%     figure;
%     subplot(221);
%     hold on;
%     plot((1:length(totSum))./length(totSum), totSum, 'LineWidth',2);
%     plot([0 t50]./length(totSum), [0.5 0.5], 'k', 'LineWidth', 2);
%     plot([t50 t50]./length(totSum), [0 0.5], 'k', 'LineWidth', 2);
%     plot([0 1], [0 1], 'k--');
%     xlim([0 1]);
%     ylim([0 1]);
%     hold off;
%     subplot(222);
%     hold on;
%     plot((1:length(outSum))./length(outSum), outSum, 'LineWidth',2);
%     plot([0 o50]./length(outSum), [0.5 0.5], 'k', 'LineWidth', 2);
%     plot([o50 o50]./length(outSum), [0 0.5], 'k', 'LineWidth', 2);
%     plot([0 1], [0 1], 'k--');
%     xlim([0 1]);
%     ylim([0 1]);
%     hold off;
%     subplot(223);
%     hold on;
%     plot((1:length(inSum))./length(inSum), inSum, 'LineWidth',2);
%     plot([0 i50]./length(inSum), [0.5 0.5], 'k', 'LineWidth', 2);
%     plot([i50 i50]./length(inSum), [0 0.5], 'k', 'LineWidth', 2);
%     plot([0 1], [0 1], 'k--');
%     xlim([0 1]);
%     ylim([0 1]);
%     hold off;
%     subplot(224);
%     hold on;
%     plot((1:length(allsyns))./length(allsyns), allsyns, 'LineWidth',2);
%     plot([0 a50]./length(allsyns), [0.5 0.5], 'k', 'LineWidth', 2);
%     plot([a50 a50]./length(allsyns), [0 0.5], 'k', 'LineWidth', 2);
%     plot([0 1], [0 1], 'k--');
%     xlim([0 1]);
%     ylim([0 1]);
%     hold off;
end

