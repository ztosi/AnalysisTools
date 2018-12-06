mnSyn11_51 = zeros(sum(inset11{51}),1);
totSyn11_51 = zeros(sum(inset11{51}),1);
mnUnq11_51 = zeros(sum(inset11{51}),1);
totUnq11_51 = zeros(sum(inset11{51}),1);
mnRed11_51 = zeros(sum(inset11{51}),1);
totRed11_51 = zeros(sum(inset11{51}),1);
iset = find(inset11{51});
for ii=1:length(iset)
[mnSyn11_51(ii), totSyn11_51(ii), mnUnq11_51(ii), totUnq11_51(ii), mnRed11_51(ii), totRed11_51(ii)] = findSynVals(51, iset(ii), pilData{11});
end