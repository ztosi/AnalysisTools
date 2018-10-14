function [symbols, freq] = spikeCompressionBase(asdf)
  
freq = cellfun(@length, asdf(1:end-2));
freq = freq ./ sum(freq);

symbols = java.huffman(freq);

%symbols = javaMethod('huffman', 'Huffman', freq);
    
    



end