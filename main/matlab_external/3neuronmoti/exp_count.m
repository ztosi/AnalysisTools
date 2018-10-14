function [N_exp] = exp_count(N,p_nc,p_uni,p_bi,song_order)
% calculate analytical expressions for 3 neuron motifs based on figure 4A 
% in Song et al.
% Input N = # of nodes. p_nc = probablity of un-connected pair evaluated 
% from actual data. p_uni = probablity of uni-directional pair; 
% p_bi = probablity of bi-directional pair. song_order = order in which 
% motifs are listed in Song's paper. 
% OUTPUT : N_exp = frequency of occurence of 13 motifs as numbered in
% Song's paper

% the numbering of motifs used to calculate probabilities is 
% initially chosen as defined in Sporns & Kotter PLoS Biology 2004 "Motifs in Brain Networks

p(1) = 3*(p_uni^2)*p_nc; p(2) = 2*p(1); p(3) = p(1);
p(4) = 6*(p_uni)*(p_bi)*(p_nc); p(5) = 6*(p_uni^3);
p(6) = p(4);
p(7) = 2*(p_uni^3); p(8)= 3*(p_uni^2)*(p_bi);
p(9) = 3*((p_bi)^2)*(p_nc); p(10)= 2*p(8); p(11) = p(8);
p(12) = 6*(p_bi^2)*(p_uni);
p(13) = p_bi^3;
% renumber to match up with Songs plot
p = p(song_order);
N_exp = (N*(N-1)*(N-2)).*p;
end

