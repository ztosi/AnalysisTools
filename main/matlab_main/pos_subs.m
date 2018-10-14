function [ s, okvals ] = pos_subs(sz, y, x )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    nn = x>0 & y>0 & x<=sz(2) & y<=sz(1);
    x=x(nn);
    y=y(nn);
    s=sub2ind(sz, y, x)';
    okvals=nn;

end

