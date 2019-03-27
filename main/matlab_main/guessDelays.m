function [eds, dlyData] = guessDelays(wtMat,  asdf,  ranges, varargin)

narginchk(3,11);

method = 'mean'; % How to determine the delay
%	pdM = 'none'; % How to determine the timing distribution
metaFlag = true;
binUnit = 1;

for vi = 1:2:length(varargin)
    switch varargin{vi}
        case 'Method'
            method = varargin{vi+1};
            %			case 'PdMethod'
            %				pdM = varargin{vi+1};
        case 'Meta'
            metaFlag = varargin{vi+1};
        case 'BinUnit'
            binUnit = varargin{vi+1};
        otherwise
            error('Unrecognized input');
    end
end

fun = @mean;
switch method
    case 'mean'
        fun = @mean;
    case 'median'
        fun = @median;
    case 'mode'
        fun = @mode;
    otherwise
        error('Unrecognized method');
end

if metaFlag
    rast = ASDFToRaster(asdf, 'BinUnit', binUnit, 'col');
else
    rast = ASDFToRaster(asdf, 'BinUnit', binUnit, 'col', '-nometa');
end

numDlys = (ranges(2)-ranges(1))./binUnit;
numDlys = ceil(numDlys);
ranges = ceil(ranges./binUnit) .* binUnit;
initDly = ranges(1)/binUnit;

inpData = cell(size(wtMat,2),1);
indData = cell(size(wtMat,2),1);

for ii=1:size(wtMat,2)
    disp(ii);
    tSpks = find(rast(:, ii));
    
    inpSet = wtMat(:,ii)~=0;
    inpRast = ASDFToRaster(asdf(inpSet), 'BinUnit', binUnit, 'row', '-nometa');
    if isempty(inpRast)
        continue;
    end
    ed = size(inpRast,2);
    
    dlyVals = zeros(nnz(inpSet), numDlys);
    
    targIsis = log(diff(tSpks) .* binUnit);
    mn = min(targIsis);
    mx = max(targIsis);
    targPd = fitdist(targIsis, 'kernel');
    band = targPd.BandWidth;
    targPdf = pdf(targPd, (mn-band/2):(mx-mn + band)/1000:(mx+band/2));
    targPdf = targPdf ./ sum(targPdf);
    targEnt = -sum(nonzeros(targPdf).*log2(nonzeros(targPdf)));
    
    tSpks = tSpks(tSpks>ranges(2));
    tSpks = tSpks(tSpks < ed);
    
    for kk=1:numDlys
        if nnz(inpSet) == 0
            break;
        end
        shift  = initDly + kk - 1;
        counts = sum(inpRast(:, tSpks-shift),2);
%         inpRast(:, tSpks-shift) = 0;
        dlyVals(:, kk) = counts;
    end
    
    allDat = zeros(nnz(inpSet),4);
    
    for jj=1:nnz(inpSet)
        data = zeros(sum(dlyVals(jj,:)), 1);
        cVal = 1;
        for kk=1:numDlys
            data(cVal:(cVal + dlyVals(jj,kk)-1)) = (kk+initDly-1) * binUnit;
            cVal = cVal + dlyVals(jj,kk);
        end
        
        if isempty(data)
            disp('data was empty');
            continue;
        end
        
        pd = fitdist(data, 'kernel');
        bw = pd.BandWidth;
        dst = pdf(pd, (ranges(1)-bw/2):(diff(ranges)+bw)/1000:(ranges(2)+bw/2));
        dst = dst ./ sum(dst);
        
        allDat(jj,1) = fun(data);
        allDat(jj,2) = -sum(nonzeros(dst).*log2(nonzeros(dst)));
        allDat(jj,3) = targEnt;
    end
    
    allDat(:,4) = allDat(:,3) - allDat(:,2);
    
    inpData{ii} = allDat;
    indData{ii} = [find(inpSet), ones(nnz(inpSet), 1) .* ii];
    
end

eds = cell2mat(indData);
dlyData = cell2mat(inpData);

end