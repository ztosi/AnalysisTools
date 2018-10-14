function [ crast ] = quickSpkTrainConv( rast, compress, tc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[m,n]=size(rast);

crast = zeros(m, ceil(n/compress));
%rast=rast';
clength = 5*tc;
kernel = exp(-(0:clength)/tc);

    for i=1:m
       if (issparse(rast))
           rw = conv(full(rast(:,i))', kernel, 'full');
       else
           rw = conv(rast(i,:), kernel, 'full');         
       end
       crast(i,:) = rw(1:compress:n);  
    end


end

