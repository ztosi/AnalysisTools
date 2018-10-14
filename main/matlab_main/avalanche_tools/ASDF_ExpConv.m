function [raster, binunit] = ASDF_ExpConv(asdf, varargin)
%tic;
narginchk(1, 3);
if ~isempty(asdf{end})
    info = asdf{end};
    n_neu = info(1);
    binunit = asdf{end - 1};
else
    n_neu = length(asdf)-2;
    binunit = 0.25;
end
% change to column vectors if need be
asdf = cellfun(@(x) reshape(x, length(x), 1), ...
    asdf, 'UniformOutput', false);
colMajor = true;
if isempty(varargin)
    binsize = binunit;
else
    for ii=1:length(varargin)
        if isnumeric(varargin{ii})
            binsize = varargin{ii};
            if binsize < binunit
                error('Cannot have custom bin size less than asdf time bin.');
            end
        else
            if strcmp(varargin{ii}, 'row') || ...
                    strcmp(varargin{ii}, 'row major')
                colMajor = false;
            elseif strcmp(varargin{ii}, 'column') || ...
                    strcmp(varargin{ii}, 'col')  || ...
                    strcmp(varargin{ii}, 'column major')
                colMajor = true;
            else
                error('Unknown input');
            end
        end
    end
end
nemp = find(cellfun(@isempty, asdf(1:end-2)));
nemp = [setdiff(1:n_neu, nemp) n_neu+1 n_neu+2];
asdf = asdf(nemp);
n_neu = length(nemp)-2;

minVal = min(cellfun(@min, asdf(1:end-2)));

% very simple check of validity of ASDF
if n_neu ~= size(asdf,1) - 2
    error('Invalid n_neu information is contained in this ASDF');
end

% Get rid of unnecessary metadata
asdf = asdf(1:(length(asdf)-2));

asdf_exp = cell(923,1);
for ii=1:length(asdf)
    asdf{ii} = asdf{ii}-minVal;
    asdf_exp{ii} = (binsize .* ceil(asdf{ii}/binsize)) - asdf{ii};
    asdf_exp{ii} = exp(-asdf_exp{ii}/(binsize/3));
end

convToBinInd = @(x) 1+floor((x ./ binsize));
asdf = cellfun(convToBinInd, asdf, 'UniformOutput', 0);

bins = cellfun(@unique, asdf, 'UniformOutput', 0);
vals = cellfun(@(x,y) nonzeros(accumarray(x, y)), asdf, asdf_exp, 'UniformOutput', 0);
V = cell2mat(vals);
I = cell2mat(bins);
J = cell2mat(cellfun(@(x,y) ones(length(x),1).*y, bins, ...
    num2cell((1:n_neu)'), 'UniformOutput', 0 ));

size(V)
size(I)
size(J)

if colMajor
    raster = sparse(I, J, V);
else
    raster = sparse(J, I, V);
end

end