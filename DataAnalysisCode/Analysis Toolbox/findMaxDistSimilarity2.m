function [ js_diverge ] = findMaxDistSimilarity2( dat1, dat2, smoothKernel)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if isrow(dat1)
    dat1=dat1';
end
if isrow(dat2)
    dat2=dat2';
end

dat1 = dat1(dat1<20000);
dat2 = dat2(dat2<20000);
dat1 = (dat1-mean(dat1))./std(dat1);
dat2 = (dat2-mean(dat2))./std(dat2);

% Use default matlab binning method to find bin-size that covers whole
% dataset (dat 1 & 2); TODO: Allow user to set bin method option or
% specify their own bins.
[~, ed] = histcounts([dat1; dat2]);
d1raw = histcounts(dat1, ed, 'Normalization', 'probability');
d2raw = histcounts(dat2, ed, 'Normalization', 'probability');
binsize = ed(2)-ed(1);
% Smooth distributions. TODO: allow user to not smooth in not dumb way.
d1 = conv(d1raw, smoothKernel, 'same');
d2 = conv(d2raw, smoothKernel, 'same');

clear d1raw;
clear d2raw;

% Renormalize
d1 = d1./sum(d1);%.*binsize);
d2 = d2./sum(d2);%.*binsize);

d1=d1';
d2=d2';

% For now just look at shift... try scale later...
[crossC, mxI] = max(conv(d1, flipud(d2), 'full'));

% how much to shift d1
shift = mxI - length(d1);

vec1 = zeros(3*length(d1)-1, 1);
vec2 = zeros(3*length(d1)-1, 1);

inds1 = (1:length(d1)) + length(d1)-1;
inds2 = inds1 + shift;

vec1(inds1) = d1;
vec2(inds2) = d2;

%figure; plot(vec1); hold on; plot(vec2); hold off;

p_pr = vec1;
q_pr = vec2;
m_pr = (p_pr + q_pr)./2;

pnz = p_pr ~= 0;
qnz = q_pr ~= 0;
d1 = 0.5*sum(p_pr(pnz) .* log2(p_pr(pnz)./ m_pr(pnz)));
d2 = 0.5*sum(q_pr(qnz) .* log2(q_pr(qnz) ./ m_pr(qnz)));

js_diverge = d1+d2;

end

