% Based on Sunny Nigam's Avalanche_detection.m
% 
% Takes 1ms-binned asdf data
% Produces shapes data with 2nd column in absolute time
%
% Rashid Williams-Garcia 12/11/15

function [aShapes,aLengths,aSizes] = asdf2Shapes(asdf)
    asdf = asdfCast(asdf);
    rast = ASDFToSparse(asdf);
    pop_fir = sum(rast);    %sum columns of the raster
    ind_fr = find(pop_fir==0)'; %find times of inactivity
    
    av = 1;
    for k = 1:length(ind_fr)-1
        if ind_fr(k+1)-ind_fr(k)>1  % keep length 1 avalanches
            ava = ind_fr(k)+1:ind_fr(k+1)-1;
            temp = rast(:,ava); % isolate those bins in the time raster
            [neuronIDs,Times] = find(temp~=0); % find out which neurons spike and when
            aSizes(av,1) = length(neuronIDs);
            aLengths(av,1) = length(unique(Times));
            aShapes{av,1} = horzcat(neuronIDs,reshape(ava(Times),size(neuronIDs)));% store this in a cell array
            av = av+1;
        end
    end
end