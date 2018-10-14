function [ div ] = jensen_shannon( p_dat, q_dat, varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    narginchk(2, 3);
    default_bin = 'fd';
    
    %p_dat = sort(p_dat, 'ascend');
    %q_dat = sort(q_dat, 'ascend');
    
    
    if isempty(varargin)
       [m_pr, ed, bin_m] = histcounts([p_dat q_dat], 'BinMethod', ...
           default_bin, 'Normalization', 'pdf');
    else
       [m_pr, ed, bin_m] = histcounts([p_dat q_dat], 'BinMethod', ...
           varargin{1}, 'Normalization', 'pdf');
    end
    
    [p_pr, ~, bin_p] = histcounts(p_dat, ed, 'Normalization', ...
           'pdf');
    [q_pr, ~, bin_q] = histcounts(q_dat, ed, 'Normalization', ...
           'pdf');

     m_pr = (p_pr + q_pr)./2;   
       
     pnz = p_pr ~= 0;
     qnz = q_pr ~= 0;
     d1 = 0.5*sum(p_pr(pnz) .* log2(p_pr(pnz)./ m_pr(pnz)));
     d2 = 0.5*sum(q_pr(qnz) .* log2(q_pr(qnz) ./ m_pr(qnz)));

     div = d1+d2;

end

