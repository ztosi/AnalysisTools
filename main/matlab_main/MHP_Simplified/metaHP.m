function prefFRs = metaHP(adjMat, inD, prefFRs, estFRs, eta, time)
    global N alpha lowFR beta dt;
    % Adj is supposed to be transposed

    max = 100;
    min = 0.01;
%     function mdnz = medianNz(m)
%        mdnz = zeros(size(m,2),1);
%         for ii=1:size(m,2)
%            mdnz(ii) = median(nonzeros(m(:,ii)));
%         end
%     end
%     
%     in = adjMat .* log(estFRs);
%     mn = sum(in) ./ sum(adjMat);
    
    b = 0.25; % b .* exp(-((log(prefFRs)-1)./.5).^2);
    
    mat = adjMat .* (log(estFRs'+0.00001)-log(estFRs+0.0001));
%    mat(mat(:)==0) = NaN;
    mnDf = median(mat, 'omitnan')';%sum(mat)' ./ inD;
%    mat(isnan(mat(:))) = 0;
    %  b = 0.4.*(1-exp(-abs(mnDf./2)))+0.05;%b.*(log(abs(mnDf-log(prefFRs)) + 1 ))+0.01;
    mat = 2./(b'.*b') .* mat .* exp(-((mat)./(b')).^2); %.* -tanh(adjMat.*(abs(mat)-1));
    
    %     pfg1=prefFRs.*((prefFRs>1).* -1/(max-1) + (1+1/max));
    %     pfl1 = prefFRs.*((prefFRs<=1) .* (1+min) - min);
    %mnDf = -2 .* mnDf .* exp(-(mnDf).^2) .* (inD./10);
    mnDf =  (2./(1+exp(2.*mnDf)) - 1);
    mnDf(mnDf < 0) = mnDf(mnDf<0) - 0.1;
    mnDf(mnDf > 0) = mnDf(mnDf>0) + 0.1;
    mnDf = mnDf .* (inD./10);
    dpfr = ((sum(mat, 'omitnan')'+mnDf)./inD);% - 0.01.*(log(prefFRs)-5).^-2 + 0.01.*(log(prefFRs)+5).^-2; %.* sqrt(1-(log(prefFRs)/4.5).^2);
    if mod(time, 200) == 0
        subplot(121);
        No = 213;
        title(['In Distribution and Smoth Function Vals for Neuron' num2str(No)]);
        divKernF = @(x,b) 2./(b.*b) .* x .* exp(-((x)./(b)).^2);
        inEfs = log(estFRs(adjMat(:,No) == 1));
        
        dom = -3:0.01:5;
        inEfInd = ceil(inEfs ./ 0.01) .* 0.01;
        bh = histcounts(inEfInd, [dom dom(end)]+0.005, 'Normalization', 'count');
        hold on;
        plot([dom(1) dom(end)], [ 0 0], '-k', 'LineWidth', 2);
        plot([inEfs, inEfs]', 0.1 * [zeros(1, length(inEfs)); ones(1, length(inEfs))], 'k');
        plot([log(prefFRs(No)) log(prefFRs(No))], [-.4 .4], 'LineWidth', 2);
        plot([log(estFRs(No)) log(estFRs(No))], [-.4 .4], 'Color', [.05 .55 .17], 'LineWidth', 2);    
        mnF = median(inEfs);
        plot([mnF, mnF], [-.4 .4], '--k', 'LineWidth', 3);
        tmp = (2./(1+exp(2.*(dom - mnF))) - 1);
        tmp(tmp<0)=tmp(tmp<0) - 0.1;
        tmp(tmp>0)=tmp(tmp>0) + 0.1;
        plot(dom, conv(bh, divKernF(-2:0.01:2, b), 'same' )./inD(No));
                der = (conv(bh, (divKernF(-2:0.01:2, b)),  'same')...
            +  (inD(No)./10) .* tmp)./inD(No);
        plot(dom, der, 'LineWidth', 1);
        plot(dom, conv(bh, exp(-((-2:0.01:2)/b).^2)./inD(No), 'same'), 'LineWidth', 3);
%         pt = (divKernF(log(estFRs(No)), .3) + (inD(No)./10) .* (2./(1+exp(log(estFRs(No)) - mnF)) - 1))./inD(No);
%         scatter(log(estFRs(No)), pt, 40, [0 0 0], 'filled');        
        xlabel(['Ln(firing rate), Mean: ' num2str(exp(mnF)) ' Hz, TFR: ' num2str(prefFRs(No)) ' Hz, b:' num2str(b) ]  );
        xlim([-3,5]);
        ylim([-.5, .5]);
        drawnow;
        hold off;
    end
    
    prefFRs = exp(log(prefFRs) + dt * eta * dpfr);
    
    
  
%     dfs = estFRs-estFRs';
%     
%     dfs = dfs .* adjMat;
%     
%     ltVals = (exp(-abs(((dfs>0).*dfs)./prefFRs')));
%     ltVals = sum(ltVals .*adjMat)';
%     gtVals = (exp(-abs(((dfs<0).*dfs)./prefFRs')));
%     gtVals = sum(gtVals .* adjMat)';
% 
%     f_plus = exp(-prefFRs./(beta*lowFR));
%     lf = prefFRs <= lowFR;
%     f_minus_L = prefFRs(lf)./lowFR;
%     f_minus = (1 + log(1+alpha*(prefFRs(~lf)./lowFR -1))./alpha);
%     gtVals(lf) = gtVals(lf) .* f_minus_L;
%     gtVals(~lf) = gtVals(~lf) .* f_minus;
%     ltVals = ltVals .* f_plus;
%     
%     prefFRs = prefFRs + eta * dt .* (1.5.*ltVals - gtVals) ./ inD .* ((0.25.*randn(N,1))+1);
    

end