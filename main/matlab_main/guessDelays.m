function [dlys, ent, tent, mi] = guessDelays(edges,  asdf,  ranges, varargin)

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

	dlys = zeros(size(edges));
	ent = zeros(size(edges,1),1);
	tent = zeros(size(edges,1),1);
	mi = zeros(size(edges,1),1);


	if metaFlag
		rast = ASDFToRaster(asdf, 'BinUnit', binUnit, 'col');
	else
		rast = ASDFToRaster(asdf, 'BinUnit', binUnit, 'col', '-nometa');
	end


	% Ensure things are sorted by target
	[~, tOrd] = sort(edges(:,2));
	edges = edges(tOrd,:);
	if ~isempty(ranges)
		ranges = ranges(:, tOrd);
	end
	targs = unique(edges(:,2));

	%TODO: Add reliability
	edgeNo = 1;
	for jj = 1:length(targs)
		tind = targs(jj);
		srcs = sort(edges(edges(:,2) == tind, 1));
		tSpks = find(rast(:, tind));

		targEnt = log(diff(asdf{tind}));
		mn = min(targEnt);
		mx = max(targEnt);
		targEnt = fitdist(targEnt, 'kernel');
		targEnt = pdf(mn:(mx-mn)/2000:mx, targEnt);
		targEnt = targEnt ./ sum(targEnt);
		targEnt = -sum(targEnt .* log2(targEnt));

		if ~isempty(ranges)
			locRanges = ceil(ranges(edges(:,2) == tind)./binUnit).*binUnit;

			for ii = 1:length(srcs)
				sSpks = find(rast(:, srcs(ii)));
				timeData = zeros(size(sSpks),1);
				timeDataInd = 1;
				whoswho = zeros(length(tSpks)+length(sSpks),1);
				whoswho(1:length(sSpks)) = 1;
				[~,ord] = sort([sSpks; tSpks], 'ascend');
				whoswho = whoswho(ord);

				n2check = (locRanges(2) - locRanges(1))/binUnit + 1;
				kernel = zeros(n2check + locRanges(1)/binUnit,1);
				data = [];
				%timeBins = zeros(n2check,1);
				for kk=1:n2check
					kernel(kk + locRanges(1)/binUnit) = 1; % Advance the point we want to check
					tmp = conv(whoswho, kernel, 'full'); % will be "1" if a source was offset from a target by kk+locRanges(1)/binUnit 
					data  = [data; nonzeros(tmp(1:length(whoswho)) .* (whoswho == 0)) .* ((kk-1)/binUnit + locRanges(1))];
					kernel(kk + locRanges(1)/binUnit) = 0; % cover up the old one
				end



				switch method
					case 'mean'
						dly(edgeNo) = mean(data); 
					case 'median'
						dly(edgeNo) = median(data);
					case 'mode'
						dly(edgeNo) = mode(data);

					otherwise
						error('Unrecognized method');

				end
				timeDist = fitdist(log(data), 'kernel');

				td = pdf(log(locRanges(1)):(log(locRanges(2) - locRanges(1))/2000):log(locRanges(2)), timeDist);

				td = td ./ sum(td);

				ent(edgeNo) = -sum(td .* log2(td));
				tent(edgeNo) = targEnt;

				mi = targEnt - ent(edgeNo);

				edgeNo = edgeNo + 1;

			end

			
		end

	end

end