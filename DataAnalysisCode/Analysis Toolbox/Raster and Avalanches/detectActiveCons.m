function [svid] = detectActiveCons(str_ava, adj_mat, delays, frs, lenRange, maxFrames, ts, tau, kappa)

	lens = cellfun(@(x) size(x,2), str_ava);
	[~, lenOrd] = sort(lens, 'descend');
	% TODO maybe have this key as an output arg so people can match?

	str_ava = str_ava(lenOrd);
	avinrange = cellfun(@(x) size(x,2) <= lenRange(2) && size(x,2) >= lenRange(1), str_ava);

	str_ava = str_ava(avinrange);


	disp([num2str(length(str_ava)) ' Avalanches fit this range criterion.']);

	delays = delays .* adj_mat;

	svid = zeros(maxFrames, nnz(adj_mat));

	% this (next 3 lines) is so we can figure out what index of svid we should use
	% from the coordinates of a connection in their matrix... I'm sure matlab has
	% some O(log n) lookup strategy implemented that'll be better than what I'd do
	[src, targ] = find(adj_mat);
	lookup = 1:nnz(adj_mat);
	lookup = sparse(src, targ, lookup); 


	alphaFun = zeros(uint32(3*tau/ts), 1);
	for ii=1:length(alphaFun)
		t=ii*ts;
		alphaFun(ii) = (t/tau)^kappa * exp(-abs(((t-tau)/tau)^kappa)); 
	end

	for kk=1:length(str_ava)
		ava = str_ava{kk};
		nUnq = size(ava,1);
		[mems,~, ~] = find(ava{kk});
		localSd= delays(mems, mems);
		alphaConvs = zeros(nUnq, maxFrames);

		for ii=1:nUnq
			locOuts = find(localSd(ii,:)~=0);
			if isempty(locOuts)
				continue;
			end
			alphaConvs(ii,:) = conv( full(str_ava{kk}(ii,:)), alphaFun, 'same');

			for jj=1:length(locOuts)
				outT = find(ava(mems(locOuts(jj)),:));
				outT = outT + localSd(ii, locOuts(jj));
				outT = outT(outT <= maxFrames);
				index = lookup(mems(ii), mems(locOuts(jj)));
				svid(index, outT) = svid(index, outT) + alphaConvs(ii,ouT); 
			end

		end


	end


end