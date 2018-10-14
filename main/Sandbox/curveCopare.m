function [js_diverge, raw] = curveCopare(d1, d2)

%offset = min([d1 d2]);

[~, mxLag] = max(conv(d1, flipud(d2), 'full'));

% how much to shift d1
shift = mxLag - length(d1);

% Padded vectors for comparing shifted dist to non shifted dist,
% represent probability distributions p and q to be compared
% respectivly
p_pr = zeros(3*length(d1)-1, 1);
q_pr = zeros(3*length(d1)-1, 1);

% Indices in padded vector for d1 (unshifted) and d2 (shfited)
% respectivly
inds1 = (1:length(d1)) + length(d1)-1;
inds2 = inds1 + shift;

% Populate padded vectors.
p_pr(inds1) = d1;
q_pr(inds2) = d2;

raw = sum((p_pr-q_pr).^2);

p_pr = p_pr ./ sum(p_pr);
q_pr = q_pr ./ sum(q_pr);

%figure; plot(p_pr); hold on; plot(q_pr); hold off;

% Compute Jensen-Shannon Divergence in base-2
m_pr = (p_pr + q_pr)./2;
pnz = p_pr ~= 0;
qnz = q_pr ~= 0;
p_pr = 0.5 * sum(p_pr(pnz) .* log2(p_pr(pnz) ./ m_pr(pnz)));
q_pr = 0.5 * sum(q_pr(qnz) .* log2(q_pr(qnz) ./ m_pr(qnz)));
js_diverge = p_pr + q_pr;