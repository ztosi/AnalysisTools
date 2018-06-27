function [svid] = detectActiveCons(str_ava, adj_mat, delays, lenRange, maxFrames, ts, tau)

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


kernelFun = exp(-(0:ts:ceil(4*(tau/ts))));

%alphaFun = zeros(uint32(3*tau/ts), 1);
% for ii=1:length(alphaFun)
%     t=ii*ts;
%     %alphaFun(ii) = (max([0, (t/tau)]) * exp((tau-t)/tau))^(kappa^2);
% end



for kk=1:length(str_ava)
    ava = str_ava{kk};
    alen = size(ava,2);
    if alen > maxFrames
        alen = maxFrames;
    end
    [mems,~, ~] = find(ava);
    mems = unique(mems);
    nUnq = length(mems);
    localSd = delays(mems, mems);
    alphaConvs = zeros(maxFrames, nUnq);
    
    for ii=1:nUnq
        locOuts = find(localSd(ii,:)~=0);
        if isempty(locOuts)
            continue;
        end
        %   size(full(ava(mems(ii),:)))
        %      aCon = conv( full(ava(mems(ii),1:alen)), alphaFun, 'full');
        %      alphaConvs(1:alen,ii) = aCon((ceil(tau)):(alen+ceil(tau)-1));
        aCon = conv( full(ava(mems(ii),1:alen)), kernelFun, 'full');
        alphaConvs(1:alen,ii) = aCon(1:alen);
        for jj=1:length(locOuts)
            outT = find(ava(mems(locOuts(jj)),1:alen));
            outT = outT - localSd(ii, locOuts(jj));
            outT = outT(outT <= alen);
            outT = outT(outT>0);
            index = lookup(mems(ii), mems(locOuts(jj)));
            svid(outT, index) = svid(outT, index) + alphaConvs(outT, ii);
%             try
%                 svid(outT, index) = svid(outT, index) + alphaConvs(outT, ii);
%             catch ME
%                 disp('poop');
%             end
        end
        
    end
    
    
end


end