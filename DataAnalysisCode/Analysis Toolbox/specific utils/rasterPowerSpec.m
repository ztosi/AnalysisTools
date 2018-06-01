function [pwMat] = rasterPowerSpec( signal, binsize )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    binsize = binsize/1000;
    
    pwMat = zeros(800, length(signal));
    freqs2check = 0.1:0.1:80;
    parfor i=1:800
       pwMat(i,:) = conv(signal, morlet(freqs2check(i), binsize), 'same'); 
    end

    sv = std(nonzeros(pwMat));
    figure; set(gca, 'Ydir', 'reverse'); imagesc((1:length(signal))/1000/binsize, freqs2check, ...
        pwMat, [-sv, 2*sv]);
  
    

end

