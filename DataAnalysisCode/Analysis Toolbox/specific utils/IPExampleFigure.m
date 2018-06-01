noi = 909
ef = EstFiringRates(noi);
f1=figure;
f2=figure;
set(f1, 'NextPlot', 'replace');
set(f1, 'Position', [50 50 750 540]);
%set(f2, 'NextPlot', 'replace');
%set(f2, 'Position', [50 50 750 1080]);
j=0;
incEfs = PrefFRs(wtMat(:,noi)~=0);
%incEfs = EstFiringRates(wtMat(:,noi)~=0);
allsumVals = cell(596,1);
pf = PrefFRs(noi);
%for ef=3:0.25:3.25
    j=j+1;
    
%     if (ef > 25)
%         ef = 50-ef;
%     end
    incEfsAbv = incEfs(incEfs>ef);
    incEfsBlw =incEfs(incEfs<ef);
    
    % Make sure the ranges will line up
    %ef = ceil(ef/0.01)*0.01;
    pf = ceil(pf/0.01)*0.01;
    incEfsAbv = ceil(incEfsAbv/0.01)*0.01;
    incEfsBlw = ceil(incEfsBlw/0.01)*0.01;
    
    curvesDwn=cell(sum(incEfs> ef), 1);
    curvesUp=cell(sum(incEfs< ef), 1);
    dwnRanges=cell(sum(incEfs> ef), 1);
    upRanges=cell(sum(incEfs< ef), 1);
    
    for i=1:length(incEfsBlw)
        range = incEfsBlw(i):0.01:(incEfsBlw(i)+(100*pf));
        upRanges{i} = range;
        curvesUp{i} = exp(-pf/(20)).*exp(-abs(range-incEfsBlw(i))/pf);
    end
    for i=1:length(incEfsAbv)
        range = (incEfsAbv(i)-(100*pf)):0.01:incEfsAbv(i);
        dwnRanges{i} = range;
        curvesDwn{i} = (1+(log(1+ (2.5*((pf/2)-1)))/2.5)).* ...
            exp(-abs(range-incEfsAbv(i))/pf);
    end
    clf(f1);
    figure(f1); hold on;
    pups = [];
    for i=1:length(incEfsAbv)
        pups(i)=plot([dwnRanges{i} incEfsAbv(i)], -[curvesDwn{i} 0], 'Color', ...
            [0, .45, .74], 'LineWidth', 2);
    end
    pdws = [];
    for i=1:length(incEfsBlw)
        pdws(i)=plot([incEfsBlw(i) upRanges{i}], [0 curvesUp{i}], 'Color', ...
            [.85, .33, .1], 'LineWidth', 2);
    end
    pfp=plot([pf pf], [-2 1.5], '-k', 'LineWidth', 2);
    efp=plot([ef ef], [-2 1.5], ...
        '--k', 'LineWidth', 2);
    xlim([-5 10*ceil(max(incEfs)/10)]);
    ylim([-2 1.5]);
    xlabel('Est. Firing Rate (Spks/s)');
    ylabel(['Individual $\frac{\partial\nu^p}{\partial t} \cdot' ...
    '\frac{1}{\eta_{ip}}$ Contributions'], 'Interpreter', 'latex');
    l=legend(gca, [pdws(1), pups(1), efp, pfp], ...
       '$f_+(\nu^p_j)exp(-|\hat{\nu}_j-\hat{\nu}_i|/\nu^p_i);\; i\in L_j(t)$',...
        '$f_-(\nu^p_j)exp(-|\hat{\nu}_j-\hat{\nu}_i|/\nu^p_i);\; i \in G_j(t)$',...
        '$\hat{\nu}_j$', '$\nu^p_j$');
    l.Interpreter = 'latex';
    
    %xlim([1 100]);
    %set(gca, 'XScale', 'log');
    hold off;
    
    drawnow;
    
    frame = getframe(f1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    if j==1
        imwrite(imind,cm,'~/IPExampleIndiv0.gif','gif', 'DelayTime',0, ...
            'Loopcount', inf);
    else
        imwrite(imind,cm,'~/IPExampleIndiv0.gif','gif', 'DelayTime', 0, ...
            'WriteMode', 'append');
    end
    
    minmaxup = zeros(length(upRanges), 2);
    minmaxdwn = zeros(length(dwnRanges), 2);
    
    minmaxup(:,1) = cellfun(@min, upRanges);
    minmaxup(:,2) = cellfun(@max, upRanges);
    minmaxdwn(:,1) = cellfun(@min, dwnRanges);
    minmaxdwn(:,2) = cellfun(@max, dwnRanges);
    

    
    totUpRange = min(minmaxup(:,1)):0.01:max(minmaxup(:,2));
    totDwnRange = min(minmaxdwn(:,1)):0.01:max(minmaxdwn(:,2));
%     totUpRange = floor(totUpRange/0.01)*0.01;
%     totDwnRange = floor(totDwnRange/0.01)*0.01;
%     minmaxup=floor(minmaxup/0.01)*0.01;
%     minmaxdwn=floor(minmaxdwn/0.01)*0.01;
    totRange = min(totDwnRange):0.01:max(totUpRange);
%     totRange = floor(totRange/0.01)*0.01;
    
    sumValUp = zeros(size(totUpRange));
    sumValDwn = zeros(size(totDwnRange));
    sumVal = zeros(size(totRange));
    
    for i=1:length(incEfsAbv)
        [~,start] = min(abs(totDwnRange-minmaxdwn(i,1)));
        %length(start:(start+length(curvesDwn{i})-1))
        %length(curvesDwn{i})
        sumValDwn(start:(start+length(curvesDwn{i})-1)) = ...
            sumValDwn(start:(start+length(curvesDwn{i})-1)) ...
            + curvesDwn{i};
    end
    for i=1:length(incEfsBlw)
        [~,start] = min(abs(totUpRange-minmaxup(i,1)));
        sumValUp(start:(start+length(curvesUp{i})-1)) = ...
            sumValUp(start:(start+length(curvesUp{i})-1)) ...
            + curvesUp{i};
    end
    
    [~, start] = min(abs(totRange-min(totUpRange)));
    sumVal(start:(start+length(totUpRange)-1)) = ...
        sumVal(start:(start+length(totUpRange)-1)) ...
        + sumValUp;
    [~, start] = min(abs(totRange-min(totDwnRange)));
    sumVal(start:(start+length(totDwnRange)-1)) = ...
        sumVal(start:(start+length(totDwnRange)-1)) ...
        - sumValDwn;
    
%    clf(f2);
    figure(f2);
    ax2=subplot(211);
    hold on;
    pup=plot(totUpRange, sumValUp, 'Color', [.85, .33, .1], 'LineWidth', 2);
    pdwn=plot(totDwnRange, sumValDwn, 'Color', [0, .45, .74], 'LineWidth', 2);
    pfp=plot([pf pf], [min([sumValUp sumValDwn])*1.2 max([sumValUp sumValDwn])*1.2], ...
        '-k', 'LineWidth', 2);
    efp=plot([ef ef], [min([sumValUp sumValDwn])*1.2 max([sumValUp sumValDwn])*1.2], ...
        '--k', 'LineWidth', 2);
    ylim([min([sumValUp sumValDwn])*1.2 max([sumValUp sumValDwn])*1.2]);
    xlim([-5 10*ceil(max(incEfs)/10)]);
    xlabel('Est. Firing Rate (Spks/s)');
    ylabel(['Sum +/- $\frac{\partial\nu_p}{\partial t} \cdot' ...
    '\frac{1}{\eta_{ip}}$ Contributions'], 'Interpreter', 'latex');
    l=legend(ax2, [pdwn, pup, efp, pfp], ...
        '$\sum_i f_+(\nu_p(i))exp(-|\nu_e(j)-\nu_e(i)|/\nu_p(i))$',...
        '$\sum_i f_-(\nu_p(i))exp(-|\nu_e(j)-\nu_e(i)|/\nu_p(i))$', '$\nu_e(i)$', '$\nu_p(i)$');
    l.Interpreter = 'latex';
    %xlim([0.1 100]);
    %set(gca, 'XScale', 'log');
    hold off;
    ax3=subplot(212);
    hold on;
    sumVal = [0, sumVal] + [sumVal, 0];
    sumVal = sumVal(1,end-1);
    sumVal = sumVal ./ 2;
    sumVal = sumVal .* ((max(totRange)-min(totRange))/length(totRange));
    sumVal = cumsum(sumVal);
    rg = max(abs(min(sumVal)*1.2), abs(max(sumVal)*1.2));
    tfp=plot(totRange, sumVal, 'Color', [.94, .69, .13], 'LineWidth', 2);
    pfp=plot([pf pf], [-rg rg], '-k', 'LineWidth', 2);
    efp=plot([ef ef], [-rg rg], ...
        '--k', 'LineWidth', 2);
    
    plot([-5 10*ceil(max(incEfs)/10)], [0 0], 'k');
    %plot([1 100], [0 0], 'k');
    xlim([-5 10*ceil(max(incEfs)/10)]);
    %xlim([1 100]);
    %set(gca, 'XScale', 'log');

    ylim([-rg, rg]);
    xlabel('Est. Firing Rate (Spks/s)');
    ylabel(['$\frac{\partial\nu_p}{\partial t} \cdot' ...
    '\frac{1}{\eta_{ip}}$'], 'Interpreter', 'latex');
    l=legend(ax3, [tfp, efp, pfp], ...
       '$\sum_{i} f_{+/-}(\nu_p(i))exp(-|\nu_e(j)-\nu_e(i)|/\nu_p(i))$',...
         '$\nu_e(i)$', '$\nu_p(i)$');
    l.Interpreter = 'latex';
    hold off;
        drawnow;
    
    frame = getframe(f2);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    if j==1
        imwrite(imind,cm,'~/IPExampleSum.gif','gif', 'DelayTime',0, ...
            'Loopcount', inf);
    else
        imwrite(imind,cm,'~/IPExampleSum.gif','gif', 'DelayTime', 0, ...
            'WriteMode', 'append');
    end

%end
