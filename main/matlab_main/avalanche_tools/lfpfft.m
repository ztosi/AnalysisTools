function [p1, f, p2, Y] = lfpfft(lfps, sampRate)

    T = 1/sampRate;
    siglen = size(lfps,1);
    t = (0:siglen).*T;
    
    
    Y = zeros(siglen, size(lfps,2));
    p2 = zeros(siglen, size(lfps,2));
    
    for ii=1:size(lfps,2)
       
        Y(:, ii) = fft(lfps(:,ii));
        p2(:,ii) = abs(Y(:,ii)./siglen);
    end
    
    p1 = p2(1:siglen/2+1,:);
    p1(2:end-1,:) = 2*p1(2:end-1,:);
    f = sampRate*(0:(siglen/2))/siglen;
    

end

