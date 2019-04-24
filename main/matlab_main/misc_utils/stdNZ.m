function [ out ] = stdNZ( mat, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    narginchk(1,2);
    if (isempty(varargin))
        dim = 1;
    else
        if isscalar(varargin{1})
            dim = varargin{1};
            if dim > ndims(mat)
                error('Requested dimension exceeds matrix dimensions');
            end
        else
            error ('Dimension input must be scalar.');
        end
    end
    mat(isnan(mat)) = 0;
    [~,J,V]  = find(mat);
    n = zeros(J(end)+1,1);
    n(2:end) = histcounts(J, J(end));
    n(1) = 1;
    n=cumsum(n);
    out = zeros(J(end),1);
    for ii = 1:J(end)        
        out(ii)=std(V(n(ii):(n(ii+1)-1)));
    end

end