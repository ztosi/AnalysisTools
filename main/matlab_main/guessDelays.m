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
		tSpks = tSpks(tSpks>ranges(2));
		inpSet = wtMat(:,ii)~=0;
		inpRast = ASDFToRaster(asdf(inpSet), 'BinUnit', binUnit, 'row', '-nometa');
		dlyVals = zeros(inpSet, numDlys);

		targIsis = diff(tSpks) .* binUnit;
		mn = min(targIsis);
		mx = max(targIsis);
		targPd = fitdist(targIsis, 'kernel');
		band = targPd.Bandwidth;
		targPdf = pdf(targPd, (mn-band):(mx-mn + 2 * band)/1000:(mx+band));
		targPdf = targPdf ./ sum(targPdf);
		targEnt = -sum(targPdf.*log2(targPdf));

		for kk=1:numDlys
			shift  = initDly + kk - 1;
			counts = sum(inpRast(:, tSpks-shift),2);
			dlyVals(:, kk) = counts;
		end

		allDat = zeros(length(inpSet),4);

		for jj=1:length(inpSet)
			data = zeros(sum(dlyVals(jj,:)), 1);
			cVal = 1;
			for kk=1:numDlys
				data(cVal:(cVal + dlyVals(jj,kk)-1)) = (kk+initDly-1) * binUnit;
			end

			pd = fitdist(data, 'kernel');
			dst = pdf(pd, 0:1.5*ranges(2)/1500:(1.5 * ranges(2)));
			dst = dst ./ sum(dst);

			allDat(jj,1) = fun(data);
			allDat(jj,2) = -sum(dst.*log2(dst));
			allDat(jj,3) = targEnt;
		end

		allDat(:,4) = allDat(:,2) - allDat(:,3);

		inpData{ii} = allDat;
		indData{ii} = [find(inpSet), ones(nnz(inpSet), 1) .* ii]; 

	end

	eds = cell2mat(indData);
	dlyData = cell2mat(inpData);
	
end