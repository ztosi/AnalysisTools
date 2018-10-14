function [ava_data] = AvaStats( str_avas, N, minAva, initBorder, dt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    initiator = zeros(N,1);
    numAppear = zeros(N,1);
    timeAppearRaw = zeros(N,1);
    timeAppearRat = zeros(N,1);
    aSizeInit = zeros(N,1);
    aSizeAll = zeros(N,1);
    aDurAll = zeros(N,1);
    aDurInit = zeros(N,1);
    inds = (1:N)';
    edges = [0;inds]+0.5;
    for ii=1:length(str_avas)
       ava = str_avas{ii};
       if size(ava,2) < minAva
          continue; 
       end
       [J, I] = find(ava');
       whoSpk = unique(I);
       jInit = J<initBorder;
       initer = unique(I(jInit));
       initiator(initer) = initiator(initer) + 1;
       aSizeInit(initer) = aSizeInit(initer) + nnz(ava);
       aSizeAll(whoSpk) = aSizeAll(whoSpk) + nnz(ava);
       aDurInit(initer) = aDurInit(initer) + size(ava,2);
       aDurAll(whoSpk) = aDurAll(whoSpk) + size(ava,2);
       numAppear = numAppear + histcounts(I, edges)'; 
       tAp = accumarray(I,J);
       tAp = tAp(whoSpk);
       timeAppearRaw(whoSpk) = timeAppearRaw(whoSpk) + tAp.*dt;
       timeAppearRat(whoSpk) = timeAppearRat(whoSpk) + tAp./size(ava,2);
    end
    timeAppearRaw = timeAppearRaw ./ numAppear;
    timeAppearRat = timeAppearRat ./ numAppear;
    aSizeInit = aSizeInit ./ initiator;
    aSizeAll = aSizeAll ./ numAppear;
    aDurInit = aDurInit ./ initiator;
    aDurAll = aDurAll ./ numAppear;
    
    stdTimeAppearRaw = zeros(N,1);
    stdTimeAppearRat = zeros(N,1);
    stdSizeInit = zeros(N,1);
    stdSizeAll = zeros(N,1);
    stdDurAll = zeros(N,1);
    stdDurInit = zeros(N,1);
    
    for ii=1:length(str_avas)
       ava = str_avas{ii};
       if size(ava,2) < minAva
          continue; 
       end
       [J, I] = find(ava');
       whoSpk = unique(I);
       jInit = J<initBorder;
       initer = unique(I(jInit));
       stdSizeInit(initer) = stdSizeInit(initer) + (nnz(ava)-aSizeInit(initer)).^2;
       stdSizeAll(whoSpk) = stdSizeAll(whoSpk) + (nnz(ava)-aSizeAll(whoSpk)).^2;
       stdDurInit(initer) = stdDurInit(initer) + (size(ava,2)-aDurInit(initer)).^2;
       stdDurAll(whoSpk) = stdDurAll(whoSpk) + (size(ava,2)-aDurAll(whoSpk)).^2;
       Jiaw = (J.*dt - timeAppearRaw(I)).^2;
       Jiat = (J./size(ava,2) - timeAppearRat(I)).^2;
       Jiaw = accumarray(I, Jiaw);
       Jiat = accumarray(I, Jiat);
       stdTimeAppearRaw(whoSpk) = stdTimeAppearRaw(whoSpk) + Jiaw(whoSpk);
       stdTimeAppearRat(whoSpk) = stdTimeAppearRat(whoSpk) + Jiat(whoSpk);
    end
    
    stdTimeAppearRaw = sqrt(stdTimeAppearRaw ./ numAppear);
    stdTimeAppearRat = sqrt(stdTimeAppearRat ./ numAppear);
    stdSizeInit = sqrt(stdSizeInit ./ initiator);
    stdSizeAll = sqrt(stdSizeAll ./ numAppear);
    stdDurInit = sqrt(stdDurInit ./ initiator);
    stdDurAll = sqrt(stdDurAll ./ numAppear);

    ava_data = struct('NumInitiated', initiator, 'NumAppear', numAppear, 'MeanTOA', ...
        timeAppearRaw, 'MeanRelTOA', timeAppearRat, 'MeanSizeOfInitiated', aSizeInit, ...
        'MeanSizeWhereAppeared', aSizeAll, 'MeanDurOfInitiated', aDurInit, ...
        'MeanDurWhereAppeared', aDurAll,  'StdTOA', stdTimeAppearRaw, ...
        'StdRelTOA', stdTimeAppearRat, 'StdSizeOfInitiated', stdSizeInit, ...
        'stdSizeWhereAppeared', stdSizeAll, 'StdDurOfInitiated', stdDurInit, ...
        'StdDurWhereAppeared', stdDurAll);
    
end

