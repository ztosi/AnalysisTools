function [avgfirst] = avgFirstFrame(str_ava, num_neu)

avgfirst = zeros(num_neu,1);
occs = zeros(num_neu,1);
for ii=1:length(str_ava)
    [I, J] = find(str_ava{ii});
    %[I, IJ] = unique(I, 'first');
    %J=J(IJ);
    ui = unique(I);
    A = accumarray(I, J);
    avgfirst(ui) = avgfirst(ui)+A;
    occs(I)=occs(I)+1;
end

avgfirst = avgfirst./occs;

end