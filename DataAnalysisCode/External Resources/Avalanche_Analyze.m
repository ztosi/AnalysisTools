function [riEvtProp, propRi, avgfirst, dur_score, initiator_score, init_sze_score] = ...
    Avalanche_Analyze(aval, durations, meanFRs, how_rich, varargin)

    durToMid = 50;
    num = sum(durations < durToMid); %& durations >= 30);

    % Find the top 1% duration avalanches 
    d_sort = sort(durations, 'descend');
    piv = d_sort(uint32(0.01 * length(d_sort)));

    avgfirst = zeros(1, length(meanFRs));
    dur_score = zeros(1, length(meanFRs));
    initiator_score = zeros(1, length(meanFRs));
    init_sze_score = zeros(1, length(meanFRs));
    numAppear = zeros(1, length(meanFRs));
    for k = 1:length(aval)
        [tmp1,  tmp2] = find(full(aval{k}));
	tmp=[tmp1, tmp2];
        if ( max(tmp(:, 1)) < 2 || max(tmp(:, 2)) > piv ) 
            continue; 
        end
        [nodes, t_ind, ~] = unique(tmp(:, 1), 'last');
        avgfirst(nodes) = avgfirst(nodes) + tmp(t_ind, 2)'; 
        dur_score(nodes) = dur_score(nodes) + ...
            (max(tmp(:, 2)) / tmp(t_ind, 2));
        initiator_score( tmp(tmp(:, 2) == 1, 1) ) =  ...
            initiator_score(tmp(tmp(:, 2) == 1, 1)) + 1;
        init_sze_score( tmp(tmp(:, 2) == 1, 1) ) =  ...
            init_sze_score(tmp(tmp(:, 2) == 1, 1)) + length(tmp(:, 1));
        numAppear(nodes) = numAppear(nodes) + 1;
    end
    avgfirst = avgfirst ./ numAppear;
    dur_score = (dur_score ./ numAppear);
    init_sze_score = init_sze_score ./ initiator_score;

    % Remove top 1% of avalanches for the remaining analyses. 
    % Heavy, power-law tail skews results, causes figures to be 
    % massive but squeeze 99% of the avalanche statistics into a corner
    % of the plot making them ureadable
    durs = durations(durations <= piv);
    aval_piv = aval(durations <= piv);
    durations = durs;
    
    dsets = length(varargin)/2;
    
    meanFRs = meanFRs ./ 1000; % convert Spk/s to Spks/ms
    
    % # of spikes from rich neurons / total # of spikes in avalanche
    % as a ratio of what would be expected given all of their
    % firing rates
    riEvtProp = zeros(dsets, length(durations));
        
    % How many different rich-club neurons participate in the avalanche
    % relative to how many one would expect given the duration of the
    % avalanche and their individual firing rates
    propRi = zeros(dsets, length(durations));
    expPartic = zeros(dsets, length(durations));
    propRi_raw = zeros(dsets, length(durations));
   
    %first50mn = zeros(dsets, 50);
    %first50stU = zeros(dsets, 50);
    %first50stL = zeros(dsets, 50);
    
    % Rich participation/"unique" rich figure. Used to replot on same 
    % subplot.
    unqF = figure;
    first50Fig = figure;
    colors = {'c', 'r', 'b', 'g', 'y', 'm'};
    
    for i = 1:2:length(varargin)
        evtF = figure;
        riPar = varargin{i};
        name = varargin{i+1};
    
        [~, ri_inds] = sort(riPar, 'descend');
        
        rich = ri_inds(1:uint32( length(ri_inds) * how_rich ) );
        num_rich = length(rich);
       % sum(meanFRs(rich) == 0)
       % min(meanFRs(rich))
       % meanFRs(rich)
        ds_ind = floor(i/2)+1;

        
        exRi = sum(num_rich .* meanFRs(rich));
        exTot = sum(length(ri_inds) .* meanFRs);
        exRatio = exRi/exTot;
        tic
        firstDTM = zeros(num, durToMid);
        fflt = repmat(1:durToMid, num, 1);
        k = 1;
        for j = 1:length(aval)
           tmp = aval{j};
           if max(tmp(:, 2)) >= durToMid %&& max(tmp(:, 2)) >= 30
              for jj = 1:durToMid
                  firstDTM(k, jj) = sum(ismember(tmp(tmp(:, 2) == jj, 1) ...
                      , rich)) / sum(meanFRs(rich));
              end
              fflt(k, :) = (fflt(k, :) <= max(tmp(:,2)));
              k = k+1;
           end
        end
        k
        size(fflt)
        sum(sum(fflt))
        fflt = fflt ~= 0;
        parfor ii = 1:length(aval_piv)
            
            mfr_loc = meanFRs;
            tmp = aval_piv{ii};
            riEvtProp(ds_ind,ii) = sum(ismember(tmp(:, 1), rich)) ...
                / length(tmp(:, 1));
            riEvtProp(ds_ind,ii) = riEvtProp(ds_ind, ii) / exRatio;
            
            % Number of rich neurons in avalanche (no repititions)
            riInvolved = length(intersect(rich, tmp(:, 1))); 
            ti = max(tmp(:, 2));
            propRi_raw(ds_ind, ii) = riInvolved / num_rich;
            expPartic(ds_ind, ii) = sum(1-((1-mfr_loc(rich)).^ti));
            propRi(ds_ind, ii) = riInvolved / expPartic(ds_ind, ii);

        end

        durs = unique(durations);
        winSz = 5;
        d = min(durs):max(durs);
        riEP_mean = zeros(1, length(d));
        riEP_median = zeros(1, length(d));
        riEP_LowErr = zeros(1, length(d));
        riEP_HighErr = zeros(1, length(d));
        riEP_Dens = cell(2, length(d));
        mx = 0;
        
        pRi_mean = zeros(1, length(d));
        pRi_raw = zeros(1, length(d));
        
        for jj = 1:length(d);
            wb = winSz;
            wu = winSz;
            if jj <= winSz
                wb = jj - 1;
            end
            if jj > length(d) - winSz
                wu = length(d) - jj; 
            end
            av_inds_of_duration = durations >= jj-wb  ...
                & durations <= jj+wu;
  
            riEP_mean(jj) = mean(riEvtProp(ds_ind, ...
                av_inds_of_duration));
            [riEP_LowErr(jj), riEP_HighErr(jj)] =  ...
                semistd(riEvtProp(ds_ind, av_inds_of_duration));
            [CdfF,CdfX] = ecdf(riEvtProp(ds_ind, av_inds_of_duration), ...
                'Function','cdf');
            BinInfo.rule = 1;
            [~,BinEdge] = internal.stats.histbins(riEvtProp(ds_ind, ...
                av_inds_of_duration),[],[],BinInfo,CdfF,CdfX);
            [y,~] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
            riEP_median(jj) = median(riEvtProp(ds_ind, ...
                av_inds_of_duration));
            
            riEP_Dens{1, jj} = BinEdge;
            riEP_Dens{2, jj} = y;
            if (max(BinEdge) > mx)
               mx = max(BinEdge); 
            end

            pRi_mean(jj) = mean(propRi(ds_ind,...
                av_inds_of_duration));
            pRi_raw(jj) = mean(propRi_raw(ds_ind, ...
                av_inds_of_duration));
        end
        
        figure(evtF); subplot(1, 2, 1);
        hold on;
        shadedErrorBar(d, riEP_mean, [riEP_HighErr; riEP_LowErr], 'b');
        plot(d, riEP_median, ':', 'LineWidth', 2, 'Color', [.2 .5 .2]);
        plot(d, ones(1, length(d)), 'k--', 'LineWidth', 2);

        title(strcat('Spikes from high ',{' '}, name, ...
            ' relative to expected'));
        hold off;
        
        ser = (0:0.01:(mx))';
        filt =  repmat(ser, 1, length(d));
        imMat = repmat(ser, 1, length(d));

       
        for ll = 1:length(d)
            %[edge, In] = sort(riEP_Dens{1, ll}, 'ascend');
            edge = riEP_Dens{1, ll};
            y = riEP_Dens{2, ll};
            %y = y(In);
            imMat(ser >= max(edge), ll) = 0;
           for mm = 2:length(edge)
               imMat(ser >= edge(mm-1) & ser <= edge(mm), ll) = y(mm-1);
           end
        end
        imMat = imMat .* (imMat ~= filt); % get rid of any orignal values missed
        disp(name);
        
        figure(first50Fig); hold on;
        %subplot(ceil(length(varargin)/4), 2, ds_ind);
        %[first50L, first50U] = semistd(firstDTM);
        [ ~, meanNZ, stdNZU, stdNZL] = statsNZ( firstDTM );
        shadedErrorBar(1:durToMid, meanNZ, [stdNZU; stdNZL], colors{ds_ind}, 1);
        plot(1:durToMid, ones(1, durToMid), 'k--');
        hold off;
        

        figure(evtF); hold on; subplot(1, 2, 2);
        imagesc(1:length(d), 0:0.01:(mx), imMat); colorbar;
        set(gca, 'YDir', 'normal');
        title(strcat('Spikes from high ',{' '}, name, ...
            ' relative to expected'));
        hold off;
        
        figure(unqF);  subplot(ceil(length(varargin)/4), 2, ds_ind);
        hold on; plot(d, pRi_mean, 'LineWidth', 2);
        plot(d, ones(1, length(d)), '--', 'LineWidth', 2, 'Color', ...
            [.5 .5 .5]);
        plot(d, pRi_raw, 'r', 'LineWidth', 2);
        plot(durs, unique(expPartic(ds_ind, :))./num_rich, 'k-', ...
            'LineWidth', 2);
        title(name);
        hold off;
        toc
    end

end
