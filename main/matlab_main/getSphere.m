function [ X, Y, Z ] = getSphere( r, npoints )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    theta = 0:(2*pi/npoints):2*pi;
    phi = ((-0.5*pi)):(pi/npoints):(pi/2);
    r_p = r * cos(phi);
    circs = repmat(r_p, npoints+1, 1);
    X = bsxfun(@times, cos(theta'), circs);
    Y = bsxfun(@times, sin(theta'), circs);
    Z = repmat(r*sin(phi), npoints+1, 1);
    
end

