 function [ js_diverge, mse ] = findMaxDistSimilarity2( dat1, dat2, varargin )
% Finds the (maximal) similarity bwetween the probability distributions 
% represented by two sets of data. Transforms the data first based on a
% selected transform. User can specify the range of the data beforehand,
% which is useful if multiple comparisons are being made between datasets
% within some larger dataset with a specific range.

if isrow(dat1)
    dat1=dat1';
end
if isrow(dat2)
    dat2=dat2';
end
tOp = 'none';
topSpec = 0;
shift = 0;
mnmxFlag = 0;
toRemove = [];
range = [0 0];
rangeSpec = 0;
for ii=1:2:length(varargin)
    
    switch varargin{ii}
        
        case 'Transform'
            toRemove = [toRemove, ii, ii+1];
            switch varargin{ii+1}
                case 'none'
                    % Do nothing
                    case 'zscore'
    %                     sd1 = std(dat1);
    %                     sd2 = std(dat2);  
                    if rangeSpec
                        warning(['zscore transforms the two disstributions' ...
                        'according to their individual parameters. No such transform will be' ...
                        'applied to the range if the range argument.' ]);
                    end
                    dat1 = (dat1 - mean(dat1));%./(sd1);
                    dat2 = (dat2 - mean(dat2));%./(sd2);
                    range = range - max([mean(dat1), mean(dat2)]);
%                     dat1 = (dat1 - mean(dat1))./(sd1 *sd2);
%                     dat2 = (dat2 - mean(dat2))./(sd1 *sd2);
%                     range = range./(sd1*sd2);
                case 'log'
                    dat1 = log(dat1);
                    dat2 = log(dat2);
                    range = log(range);
                case 'log2'
                    dat1 = log2(dat1);
                    dat2 = log2(dat2);
                    range = log2(range);
                case 'log10'
                    dat1 = log10(dat1);
                    dat2 = log10(dat2);
                    range = log10(range);
                case 'min-max'
                    mnmxFlag = 1;
                    range = [0 1];
                    if rangeSpec
                       warning('Specified range discarded due to transformation selection'); 
                    end
                case 'Box Cox'
                    %pn = 1;
                    % TODO: Get this to return the transform function so it
                    % can be applied to the range if it was sepcified. 
                    [dat1, ~] = autoBoxCox(dat1);
                    [dat2, ~] = autoBoxCox(dat2);
                otherwise
                    error('Unrecognized input');
            end
        case 'Shift'
            toRemove = [toRemove, ii, ii+1];
            shift = varargin{ii+1};
        case 'FullRange'
            if(ii~=1)
                error('FullRange, if specified must be the first optional argument');
            end
            rangeSpec = 1;
            toRemove = [toRemove, ii, ii+1];
            range = varargin{ii+1};
        otherwise
            error('Unrecognized input');
    end
    
end

shift = shift ~= 0;

if mnmxFlag
    dat1 = (dat1-min(dat1))./(max(dat1)-min(dat1));
    dat2 = (dat2-min(dat2))./(max(dat2)-min(dat2));
    mxmx=1;
    mimi=0;
else
    if isequal(range, [0 0])
    mxmx = max([max(dat1), max(dat2)]);
    mimi = min([min(dat1), min(dat2)]);
    else
        mxmx = range(2);
        mimi = range(1);
    end
end

varargin = varargin(setdiff(1:length(varargin), toRemove));

if isempty(varargin)
    pd1 = fitdist(dat1, 'kernel', 'Kernel', 'epanechnikov');
    pd2 = fitdist(dat2, 'kernel', 'Kernel', 'epanechnikov');
else
    pd1 = fitdist(dat1, varargin);
    pd2 = fitdist(dat2, varargin);
end

%    function rmse = curveDist(mu, sig)
%         pd2 = fitdist((dat2+mu)*sig
%     
%     end

step = (mxmx-mimi)/2000;
x_vals = mimi:step:mxmx;

pdf1p = pdf(pd1, x_vals);
pdf2p = pdf(pd2, x_vals);

pdf1p = pdf1p ./ sum(pdf1p);
pdf2p = pdf2p ./ sum(pdf2p);

if shift
    % Find max cross correlation and lag where it happens...
    [~, mxLag] = max(conv(pdf(pd1, x_vals), ...
        pdf(pd2, x_vals), 'full'));
    
    % how much to shift d1
    shift = uint32((length(x_vals)-mxLag)/2);
end

if shift
    pdf1 = zeros(1, length(pdf1p)+abs(shift));
    pdf2 = zeros(1, length(pdf2p)+abs(shift));
    if shift<0
        pdf2(1:length(pdf2p)) = pdf2p;
        pdf1((-shift+1):end) = pdf1p;
    elseif shift > 0
        pdf2((shift+1):end) = pdf2p;
        pdf1(1:length(pdf1p)) = pdf1p;
    else
        pdf1 = pdf1p;
        pdf2 = pdf2p;
    end
else
    pdf1 = pdf1p;
    pdf2 = pdf2p;
end

p1nz = pdf1~=0;
p2nz = pdf2~=0;
m = (pdf1+pdf2)/2;
d1 = 0.5*sum(pdf1(p1nz) .* log2(pdf1(p1nz)./m(p1nz)));
d2 = 0.5*sum(pdf2(p2nz) .* log2(pdf2(p2nz)./m(p2nz)));

js_diverge = d1+d2;

mse = mean((pdf1p-pdf2p).^2);

end

