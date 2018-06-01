f = figure;
set(f, 'Position', [50 50 1500 650]);
set(f, 'NextPlot', 'replace');
%bins = ylimD:0.25:ylimU;
interval = 2;
i=t;
    clf;
    subplot(2, 3, [1 2]);
    hold on;
    xlabel('Time (s)');
    ylabel('Target Firing Rate (Spks/s)');
   % ylim([ylimD ylimU]);
    xlim([1 t]);
    for j=1:length(FiringRates)
        if ei(j) == 1
            patchline([2:i],(EFRs(2:i, j))', 'edgecolor',...
                [0.85 0.33 0.1], 'edgealpha', 0.25);
        else
            patchline([2:i], (EFRs(2:i, j))', 'edgecolor', ...
                [0 0.43 0.75], 'edgealpha', 0.25);
        end
    end


    % ...
        %'Color', [0 0.43 0.75]);
    %plot([i i], [ylimD ylimU], 'k', 'LineWidth', 3);
    set(gca, 'YScale', 'log');
    ylim([0.01, 100]);
    %lim=get(gca, 'YLim');
    hold off;
    
    subplot(2, 3, 3);
    %hold on;
    %h=overlayedHistogram(log(PFs(i,~ei)), log(PFs(i,ei)), 20, 'Normalization', 'probability');
    h=histogram(log10(PFs(i,:)), 'Normalization', 'probability' );
    %[bh, be] = histcounts(log(PFs(i,:)), bins, 'Normalization', 'probability' );
    %sz = be(2)-be(1);
    %p=patch([be(1)-sz be(1:(end-1))+(sz/2) be(end)+sz], [0 bh 0 ], [0 0.43 0.75]);
    h.EdgeColor = 'none';
    h.FaceAlpha = 0.7;
    %for k=1:length(h)
    %    h(k).EdgeColor = 'none';
    %end
    xlim([-2 2]);
    %h.FaceAlpha = 0.75;
   % xlabel('Log(FiringRates (Hz))');
    ylabel('p(\nu_p)');
   % xlim([ylimD ylimU]);
   % ylim([0 .15]);
    view(90,-90);
%    title(['Time: ' num2str(i-2) 's; ' num2str(100*nnz(flipped)/length(flipped)) '%']);
    %hold off;
     subplot(2, 3, [4 5]);
    hold on;
    xlabel('Time (s)');
    ylabel('Threshold (mV)');
   % ylim([ylimD ylimU]);
    xlim([1 t]);
    for j=1:length(FiringRates)
        if ei(j) == 1
            patchline([2:i],(thresholds(2:i, j))', 'edgecolor',...
                [0.85 0.33 0.1], 'edgealpha', 0.25);
        else
            patchline([2:i], (thresholds(2:i, j))', 'edgecolor', ...
                [0 0.43 0.75], 'edgealpha', 0.25);
        end
    end


    % ...
        %'Color', [0 0.43 0.75]);
 %   plot([i i], [ylimD ylimU], 'k', 'LineWidth', 3);
    hold off;
    
    subplot(2,3,6);
   % ylim([ylimD ylimU]);
    %hold on;
    %h=overlayedHistogram(log(PFs(i,~ei)), log(PFs(i,ei)), 20, 'Normalization', 'probability');
    h=histogram(thresholds(i,:), 'Normalization', 'probability' );
    %[bh, be] = histcounts(log(PFs(i,:)), bins, 'Normalization', 'probability' );
    %sz = be(2)-be(1);
    %p=patch([be(1)-sz be(1:(end-1))+(sz/2) be(end)+sz], [0 bh 0 ], [0 0.43 0.75]);
    h.EdgeColor = 'none';
    h.FaceAlpha = 0.7;
    %for k=1:length(h)
    %    h(k).EdgeColor = 'none';
    %end
    
    %h.FaceAlpha = 0.75;
    xlabel('Thresholds (mV)');
    ylabel('p(\theta)');
    
   % ylim([0 .15]);
    view(90,-90);
    %title(['Time: ' num2str(i-2) 's; ' num2str(100*nnz(flipped)/length(flipped)) '%']);