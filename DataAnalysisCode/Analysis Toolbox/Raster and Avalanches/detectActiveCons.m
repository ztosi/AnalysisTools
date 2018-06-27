function [svid] = detectActiveCons(str_ava, adj_mat, delays, lenRange, maxFrames, ts, tau)
%% [svid] = detectActiveCons(str_ava, adj_mat, delays, lenRange, maxFrames, ts, tau)
%
%   [INPUTS]
%
% str_ava - a cell array of neuronal avalanches stored as sparse matrices of size N x N where 
%   N is the total number of neurons across avalanches and the same as the number of rows and
%   columns of adj_mat and where rows correspond to neurons and columns to time bins
%
% adj_mat - an adjacency matrix that is 1 where there is a connection from the neuron at row 
%   i to the neuron at column j and 0 otherwise 
%
% delays - the matrix of temporal delays of each connection 
%
% lenRange - the range of avalanche durations to run the detection on e.g. [1 20]
%   would perform this on all avalanches of duration 1 time bin to duration 20 time bin. Don't 
%   actually do this though, aside from it being dumb to run code that detects what synapses
%   are active in an avalanche that's only 1 bin unit in duration, they tend to have a 
%   power-law distribution so there's going to be a crap-ton of them... In any case. All
%   results will be aligned temporally, so that if say your longest avalanches are 20, but some
%   are 15, those avalanches will be compared matching the first to the first and so on such that
%   the last 5 frames will only give results for the ones of duration 20.
%
% maxFrames - the maximum number of time-bins to actually consider. For instance you may
%   be interested in avalanches between 150-200 time bins long, but only care what happens 
%   during the first 50 time bins
%
% ts - the time-step of the avalanches in str_ava. This is the temporal width of the
%   aforementioned "time bin"
%
% tau - detection uses an exponential kernel over the spike trains of the sending neurons
%   this allows for some wiggle room in terms of detecting if a connection got used that both
%   accounts for our own fuzziness on the actual time delay of a connection, but also that post-synaptic
%   potentials don't occur instantaneously and thus a spike that arrived at a given time might still influence
%   the post-synaptic neuron into the future, though that influence diminishes exponentially
%
%   [OUTPUTS]
%
% svid - a matrix of dimension maxFrames x M where M is the number of nonzero entries in adj_mat. 
%   Each row of svid corresponds to a "frame" a single time-bin in the avalanche while each column
%   corresponds to a synapse in a "flattened" order. If one were to do [I,J] = find(adj_mat); the 
%   columns of svid would be arranged in the same order and have the same number of elements as I or J.
%   Thus each row gives the total activity over each connection across all avalanches at a given point
%   in time. Each column gives the total activity over a given connection across all avalanches over
%   the duration of the avalanche. 
%
% [AUTHOR]
% Zoe Tosi (ztosi@iu.edu)
%
% TODOS: 6/27/18
%   1) Clean up args, use varargin w/ defaults
%   2) Allow users to select arbitrary kernels
%   3) Return the connection activity for each individual avalanche in addiion to the averaged value
%   4) comment more...
%

lens = cellfun(@(x) size(x,2), str_ava);
[~, lenOrd] = sort(lens, 'descend');
% TODO maybe have this key as an output arg so people can match?
str_ava = str_ava(lenOrd);

% Get only the avalanches with the specified durations
avinrange = cellfun(@(x) size(x,2) <= lenRange(2) && size(x,2) >= lenRange(1), str_ava);
str_ava = str_ava(avinrange);

disp([num2str(length(str_ava)) ' Avalanches fit this range criterion.']);

% Ensure the delay matrix doesn't have more nonzero entries than adj_mat
delays = delays .* adj_mat;

% Initialize the output
svid = zeros(maxFrames, nnz(adj_mat));

% this (next 3 lines) is so we can figure out what index of svid we should use
% from the coordinates of a connection in their matrix... I'm sure matlab has
% some O(log n) lookup strategy implemented that'll be better than what I'd do
% To be clear I'm using a sparse matrix to create a lookup table so I can 
% easily translate between subscript and linear indices
[src, targ] = find(adj_mat);
lookup = 1:nnz(adj_mat);
lookup = sparse(src, targ, lookup);

% create the kernel function... TODO: see TODOs
kernelFun = exp(-(0:ts:ceil(4*(tau/ts))));


for kk=1:length(str_ava)
    ava = str_ava{kk};
    alen = size(ava,2);
    if alen > maxFrames
        alen = maxFrames;
    end

    [mems,~, ~] = find(ava);
    % The neurons active in this avalanche
    mems = unique(mems);
    nUnq = length(mems);

    % A local delay matrix comprised of the sub-matrix
    % of the delay matrix made from only neurons participating
    % in this avalanche
    localSd = delays(mems, mems);

    % TODO: Wait a sec, does this need to be more than just a single column,
    % why store all the covolutions, they don't get used again...
    kConvs = zeros(maxFrames, nUnq);.
    
    for ii=1:nUnq
        locOuts = find(localSd(ii,:)~=0);
        if isempty(locOuts)
            continue;
        end
        % Convolve the sending spike train with the kernel
        aCon = conv( full(ava(mems(ii),1:alen)), kernelFun, 'full'); 
        % align it and store
        kConvs(1:alen,ii) = aCon(1:alen);

        % iterate over all the neurons in this avalanche that recieve 
        % connections from the sender ii
        for jj=1:length(locOuts)
            % Collect the spike times in this ava for the receiver
            outT = find(ava(mems(locOuts(jj)),1:alen));
            % Subtract the temporal delay of the sender->reciver connection
            % thus aligning them if the reciever spiked around when a signal
            % from the dender would've reached them
            outT = outT - localSd(ii, locOuts(jj));

            % exclude all spikes which may have ended up out of bounds
            outT = outT(outT <= alen);
            outT = outT(outT>0);

            % connection index in svid
            index = lookup(mems(ii), mems(locOuts(jj)));
            svid(outT, index) = svid(outT, index) + kConvs(outT, ii);
        end
    end   
end


end