function [ keyfr, keyno ] = frs2piano( PrefFRs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [fr, ed, keyno] = histcounts(log(PrefFRs), -1:(6/88):(5-(6/88)));
    %[fr, ed, keyno] = histcounts(PrefFRs, 1:88);
    keyfr = (2.^((keyno-49)./12)) * 440;

end

