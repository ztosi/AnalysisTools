function [kernel] = nrmFilt2D(x, y, sigx,sigy, nrm)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[X, Y] = meshgrid(x, y);

 kernel = (1/sqrt(2*sigx*pi)) .* exp(- (X.^2)/(sigx*sigx));
 kernel = kernel .* (1/sqrt(2*sigy*pi)) .* exp(- (Y.^2)/(sigy*sigy));

 if strcmp(nrm, 'normalize')
     kernel = kernel ./ sum(sum(kernel));
 elseif strcmp(nrm, 'max')
    kernel = kernel ./ max(max(kernel)); 
 end
 
end

