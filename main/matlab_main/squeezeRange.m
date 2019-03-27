function [out] = squeezeRange(x, mn, mx, datMn, datMx)

if nargin == 3
    out = (x-min(x))./(max(x)-min(x));
else
    out = (x-datMn)./(datMx-datMn);
end
out = (out.*(mx-mn)) + mn;


end