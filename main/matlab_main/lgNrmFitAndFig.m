function [ dat, vals ] = lgNrmFitAndFig(data, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    narginchk(1, 5);
    NUM_BINS = 50;
    censorLow = NaN;
    censorHigh = NaN;
    %cutoff = min(data);
    
    if nargin >=2
        NUM_BINS = varargin{1};
    end
    if nargin >=3
        censorLow = varargin{2};
    end
    if nargin >=4
        censorHigh = varargin{3};
    end
    if nargin==5
        ax = varargin{4};
    else
        figure;
        ax = gca;
    end
    hold on;
    if iscolumn(data) && ~iscell(data)
        data=data';
    end
    [m, n] = size(data);
    
    dat = [];
    if iscell(data)
        if n > 1
            error('No columnar cell arrays');
        end
        for i = 1:m
            loc = data{i};
            [r, c] = size(loc);
            if r>1
                if c>1
                    error('All data in cell arrays must be in 1D arrays.');
                end
                loc = loc';
            end
            dat = [dat, loc];
        end
    else
        dat = data(:);
    end
    dat = nonzeros(dat); % removing zeros
    dat0 = dat;
    if isnan(censorLow)
        censorLow = min(dat);
    else
        dat = dat(dat>censorLow);
    end
    if isnan(censorHigh)
        censorHigh = max(dat);
    else
        dat = dat(dat<censorHigh);
    end
    dat = dat(dat<censorHigh);
    dat = dat(dat>censorLow);
    pd = fitdist(log10(dat), 'normal')
    % pd2 = fitdist(log(dat), 'logistic')
    %bins = logspace(log10(censorLow), log10(censorHigh), NUM_BINS);
    [bh, edges] = histcounts(log10(dat), 'Normalization', 'pdf');
    %sz = edges(2)-edges(1);
    binCent = edges(1:(length(edges)-1)) + (diff(edges)./2);
    xs = linspace(min(binCent), max(binCent), NUM_BINS*100);
    dev=  std(log10(dat));
    xsplt = linspace(min(binCent)-dev, max(binCent)+dev, NUM_BINS*100);
    %histogram(dat, bins, 'Normalization','pdf');
    
    %plot(ax, xs, pdf(pd2, xs), 'LineWidth', 2);
    %bh = sz .* bh ./ diff(exp(edges));
    if m > 1
        vals = zeros(m, length(edges)-1);
        for i = 1:10
            
            if iscell(data)
                dloc = data{i};
            else
                dloc = data(i,:);
            end
            dloc = dloc(dloc<censorHigh);
            dloc = dloc(dloc>censorLow);
            [vals(i,:), ~] = histcounts(log10(dloc), edges, ...
                'Normalization', 'pdf');
            
        end
        
    else
        vals = zeros(10000, length(edges)-1);
        for i = 1:10000
            dloc = dat(rand(size(dat))<0.2);
            %dloc = dloc(dloc<censorHigh);
            %dloc = dloc(dloc>censorLow);
            [vals(i,:), ~] = histcounts(log10(dloc), edges, ...
                'Normalization', 'pdf');
        end     
    end
    
    vals(isnan(vals))=0;
    [L, U] = semistd(vals);
    axes(ax);
    errbar(binCent, mean(vals), L, U, 'k');
    mn = mean(vals);
    %scatter(ax, binCent(bh>0), mn(bh>0), 30, 'g', 'filled');
    plot(ax, xsplt, pdf(pd, xsplt), 'LineWidth', 2);
    scatter(ax, binCent(bh>0), bh(bh>0), 30, 'k', 'filled');
    %plot(ax, xs, pdf(pd2, xs), 'LineWidth', 2);
    % set(ax, 'YScale', 'log');
end

