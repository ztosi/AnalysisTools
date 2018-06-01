function [ karr, kinarr, koutarr ] = plot_degrees( adjmats, useZeros )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

karr=[];
kinarr=[];
koutarr=[];
if iscell(adjmats)
    noSets = size(adjmats,1);
    ks = cell(noSets,1);
    kIns = cell(noSets,1);
    kOuts = cell(noSets,1);
    
    for i=1:noSets
        [kIns{i}, kOuts{i}, ks{i}] = nodeDegrees(adjmats{i});
        karr = [karr ks{i}];
        kinarr = [kinarr kIns{i}];
        koutarr = [koutarr kOuts{i}];
    end
    
else
    noSets = size(adjmats,3);
    ks = cell(noSets,1);
    kIns = cell(noSets,1);
    kOuts = cell(noSets,1);
    
    for i=1:noSets
        [kIns{i}, kOuts{i}, ks{i}] = nodeDegrees(adjmats(:,:,i));
        karr = [karr ks{i}];
        kinarr = [kinarr kIns{i}];
        koutarr = [koutarr kOuts{i}];
    end
    
end

if ~useZeros
    karr = nonzeros(karr)';
    kinarr = nonzeros(kinarr)';
    koutarr = nonzeros(koutarr)';
    ks = cellfun(@nonzeros, ks, 'UniformOutput', false);
    kIns = cellfun(@nonzeros, kIns, 'UniformOutput', false);
    kOuts = cellfun(@nonzeros, kOuts, 'UniformOutput', false);
end

if noSets < 5

[kp, kb] = histcounts(karr, 'BinMethod', 'fd', 'Normalization',...
    'probability');
[kip, kib] = histcounts(kinarr, 'BinMethod', 'fd', 'Normalization',...
    'probability');
[kop, kob] = histcounts(koutarr, 'BinMethod', 'fd', 'Normalization',...
    'probability');

else
    [kp, kb] = histcounts(karr, 'BinMethod', 'integer', 'Normalization',...
    'probability');
[kip, kib] = histcounts(kinarr, 'BinMethod', 'integer', 'Normalization',...
    'probability');
[kop, kob] = histcounts(koutarr, 'BinMethod', 'integer', 'Normalization',...
    'probability');
end

if noSets > 1
    length(kp)
    kpss = zeros(noSets, length(kp));
    kipss = zeros(noSets, length(kip));
    kopss = zeros(noSets, length(kop));
    for i=1:noSets
        [kpss(i,:), ~] = histcounts(ks{i}, kb, 'Normalization', ...
            'probability');
        [kipss(i,:), ~] = histcounts(kIns{i}, kib, 'Normalization', ...
            'probability');
        [kopss(i,:), ~] = histcounts(kOuts{i}, kob, 'Normalization', ...
            'probability');
    end
    
else
    kpss = zeros(1000, length(kp));
    kipss = zeros(1000, length(kip));
    kopss = zeros(1000, length(kop));
    for i=1:1000
        [kpss(i,:), ~] = histcounts(karr(randperm(length(karr), ...
            uint32(length(karr)/10))), kb, 'Normalization', 'probability');
        [kipss(i,:), ~] = histcounts(kinarr(randperm(length(karr), ...
            uint32(length(karr)/10))), kib, 'Normalization', 'probability');
        [kopss(i,:), ~] = histcounts(koutarr(randperm(length(karr), ...
            uint32(length(karr)/10))), kob, 'Normalization', 'probability');
    end
end

kb = kb(1:end-1)+0.5;
kib = kib(1:end-1)+0.5;
kob = kob(1:end-1)+0.5;

[lowstdsK, upstdsK] = semistd(kpss);
[lowstdsKi, upstdsKi] = semistd(kipss);
[lowstdsKo, upstdsKo] = semistd(kopss);

gauss = normpdf(-10:10, 0, 2);
lowstdsK = conv(lowstdsK, gauss, 'same');
upstdsK = conv(upstdsK, gauss, 'same');
lowstdsKi = conv(lowstdsKi, gauss, 'same');
upstdsKi = conv(upstdsKi, gauss, 'same');
lowstdsKo = conv(lowstdsKo, gauss, 'same');
upstdsKo = conv(upstdsKo, gauss, 'same');
kp = conv(kp, gauss, 'same');
kp=kp./sum(kp);
kip = conv(kip, gauss, 'same');
kip=kip./sum(kip);
kop = conv(kop, gauss, 'same');
kop=kop./sum(kop);

figure; 
subplot(321);
nza = ((kp-lowstdsK) > 0) & (kp > 0) & (upstdsK>0) & (kb>0);
shadedErrorBar(kb(nza), kp(nza), [upstdsK(nza); lowstdsK(nza)], {'k-', ...
    'LineWidth', 2});
c = ceil(1+log10(max(upstdsK)));
f = floor(log10(min(lowstdsK(nza))));
ylim([10^f, 10^c]);
xlim([1 max(kb)*1.5]);
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('P(k) Log-scale');
xlabel('Degree (k) Log-scale');
subplot(322);
h=histogram(karr, kb, 'Normalization', 'probability');
h.FaceColor = [0 0 0];
h.EdgeColor = [0 0 0];
h.FaceAlpha = 0.7;
h.EdgeAlpha = 0.1;
ylabel('P(k)');
xlabel('Degree (k)');

subplot(323);
nza = ((kip-lowstdsKi) > 0) & (kip > 0) & (upstdsKi>0) & (kib>0);
shadedErrorBar(kib(nza), kip(nza), [upstdsKi(nza); lowstdsKi(nza)], {'r-', ...
    'LineWidth', 2});
c = ceil(1+log10(max(upstdsKi)));
f = floor(log10(min(lowstdsKi(nza))));
ylim([10^f, 10^c]);
xlim([1 max(kib)*1.5]);
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('P(k_{in}) Log-scale');
xlabel('In-Degree (k_{in}) Log-scale');
subplot(324);
h=histogram(kinarr, kib, 'Normalization', 'probability');
h.FaceColor = [.5 0 0];
h.EdgeColor = [.5 0 0];
h.FaceAlpha = 0.7;
h.EdgeAlpha = 0.1;
ylabel('P(k_{in})');
xlabel('Degree (k_{in})');

subplot(325);
nza = ((kop-lowstdsKo) > 0) & (kop > 0) & (upstdsKo>0) & (kob>0);
shadedErrorBar(kob(nza), kop(nza), [upstdsKo(nza); lowstdsKo(nza)], {'b-', ...
    'LineWidth', 2});
c = ceil(1+log10(max(upstdsKo)));
f = floor(log10(min(lowstdsKo(nza))));
ylim([10^f, 10^c]);
xlim([1 max(kob)*1.5]);
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('P(k_{out}) Log-scale');
xlabel('Out-Degree (k_{out}) Log-scale');
subplot(326);
h=histogram(koutarr, kob, 'Normalization', 'probability');
h.FaceColor = [0 .1 .6];
h.EdgeColor = [0 .1 .6];
h.FaceAlpha = 0.7;
h.EdgeAlpha = 0.1;
ylabel('P(k_{out})');
xlabel('Out-Degree (k_{out})');

end

