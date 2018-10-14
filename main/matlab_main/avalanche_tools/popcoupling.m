function [ popCouplings ] = popcoupling( rast, dt, kernel )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    meanFRs = sum(rast,1)./ (size(rast, 1)*dt);
    
    popVal = zeros(size(rast,1),1);
    
    for ii=1:size(rast,2)
        popVal = popVal + ...
            conv(full(rast(:,ii)), kernel, 'same')-meanFRs(ii);
    end
    popCouplings = zeros(size(rast,2),1);
    for ii=1:size(rast,2)
        popCouplings(ii) = 1/sum(rast(:,ii)).* sum(conv(full(rast(:,ii)), ...
            kernel, 'same') .* popVal);
    end

end

