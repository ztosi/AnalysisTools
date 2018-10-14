function [rastM, lfps] = showDemarcedRaster(modules, rast, cmap, varargin)
narginchk(3, 6);
[m, n]= size(rast);
ei = ones(m,1);
lineSpec = 'm-';
if ~isempty(varargin)
    ei = varargin{1};
    if length(varargin) >1 
        lineSpec = varargin{2};
    end
    if length(varargin) > 2
       axes(varargin{3}) 
    else
        figure;
    end
else
   figure; 
end
hold on;


lfps = zeros(length(modules),n);
modSzes = cellfun(@(x) length(x), modules);
    function z = makeCol(vec)
        if isrow(vec)
            vec=vec';
        end
        z = vec;
    end
modules = cellfun(@(x) makeCol(x), modules, 'UniformOutput' , 0);
modOrder = cell2mat(modules);

lineYs = cumsum(modSzes);
ei = double(ei);
ei(ei==0) = -1;
ei=makeCol(ei);

rast = rast .* ei;


colormap(cmap);
rastM = sparse(rast(modOrder,:));
imagesc(rast(modOrder,:));
lineYs = [0; lineYs];
for ii=1:(length(lineYs)-1)
   lfps(ii,:) = sum(abs(rastM((lineYs(ii)+1):lineYs(ii+1),:))); 
   plot([0 n+0.5], [lineYs(ii+1) lineYs(ii+1)], lineSpec);
end
xlim([0 n]);
ylim([0 length(modOrder)]);
hold off;

end