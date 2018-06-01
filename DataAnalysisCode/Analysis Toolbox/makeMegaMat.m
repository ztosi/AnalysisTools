function [ mat ] = makeMegaMat( inMat1, r1, ei, r12r2, r2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    m = size(r2,2) + size(r1,2) + size(inMat1,1);
    insz = size(inMat1,1);
    r1sz = size(r1,2);
    r2sz = size(r2,2);
    mat = zeros(m);
    mat(1:insz, (insz+1):(insz+r1sz)) = inMat1;
    mat((insz+1):(insz+r1sz), insz+1:(insz+r1sz)) = r1;
    ei=find(ei)+(insz); 
    mat(ei, (insz+r1sz+1):(insz+r1sz+r2sz)) = r12r2;
    mat((insz+r1sz+1):end, (insz+r1sz+1):end) = r2;
end

