function [ intraMInp, intraMOut ] = moduleFlow( wtMat, Modules )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    numNeu = length(wtMat(:,1));

    inpFrom = cell(numNeu);
    outTo = cell(numNeu);
    
    for ii=1:numNeu
       inpFrom{ii} = find(wtMat(:,ii));
       outTo{ii} = find(wtMat(ii,:));
    end
    
    intraMInp_ = zeros(numNeu, 2);
    intraMOut_ = zeros(numNeu, 2);
    
    for ii=1:length(Modules)
        for jj=1:length(Modules{ii})
           neuInd = Modules{ii}(jj);
           
           % find what proportion of the inputs to neuron neuInd come from
           % the module it belongs to
           intraMInp_(neuInd,1) = length(intersect(Modules{ii}, inpFrom{neuInd})) ...
               / length(inpFrom{neuInd});
           intraMInp_(neuInd, 2) = intraMInp_(neuInd,2) + 1;
           
           % find what proportion of the outputs from neuron neuInd go to
           % the module it belongs to
           intraMOut_(neuInd,1) = length(intersect(Modules{ii}, outTo{neuInd})) ...
               / length(outTo{neuInd});
           intraMOut_(neuInd, 2) = intraMOut_(neuInd,2) + 1;
        end
    end
    
    intraMInp_ = intraMInp_(:, 1) ./ intraMInp_(:,2);
    intraMOut_ = intraMOut_(:, 1) ./ intraMOut_(:,2);
    
    intraMInp = intraMInp_;
    intraMOut = intraMOut_;


end

