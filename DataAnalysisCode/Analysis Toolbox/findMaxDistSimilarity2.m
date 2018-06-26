function [ js_diverge ] = findMaxDistSimilarity2( dat1, dat2, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if isrow(dat1)
    dat1=dat1';
end
if isrow(dat2)
    dat2=dat2';
end
tOp = 'none';
topSpec = 0;
shift = 0;
for ii=1:length(varargin)
   if strcmp(varargin{ii}, 'Transform')
       tOp = varargin{ii+1};
       topSpec = ii+1;
   end
   
   if strcmp(varargin{ii}, 'shift')
       shift = ii;
   end
   
end

varargin = varargin(setdiff(1:length(varargin), nonzeros([topSpec topSpec-1 shift])));

shift = shift ~= 0;

%pn = 0;

mnmxFlag = 0;

switch tOp
    
    case 'none'
        % Do nothing
    case 'z-score'
        dat1 = (dat1 - mean(dat1))./std(dat1);
        dat2 = (dat2 - mean(dat2))./std(dat2);
    case 'log'
        dat1 = log(dat1);
        dat2 = log(dat2);
    case 'log2'
        dat1 = log2(dat1);
        dat2 = log2(dat2);
    case 'log10'
        dat1 = log10(dat1);
        dat2 = log10(dat2);
    case 'min-max'
        mnmxFlag = 1;
    case 'Box Cox'
        %pn = 1;
        [dat1, ~] = autoBoxCox(dat1);
        [dat2, ~] = autoBoxCox(dat2);
    otherwise
        error('Unrecognized input');
end

if mnmxFlag
    dat1 = (dat1-min(dat1))./(max(dat1)-min(dat1));
    dat2 = (dat2-min(dat2))./(max(dat2)-min(dat2));
    mxmx=1;
    mimi=0;
else
    mxmx = max([max(dat1), max(dat2)]);
    mimi = min([min(dat1), min(dat2)]);
end


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

if shift
    % Find max cross correlation and lag where it happens...
    [~, mxLag] = max(conv(pdf(pd1, x_vals), ...
        pdf(pd2, x_vals), 'full'));
    
    % how much to shift d1
    shift = uint32((length(x_vals)-mxLag)/2);
end


pdf1p = pdf(pd1, x_vals);
pdf2p = pdf(pd2, x_vals);

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

end

