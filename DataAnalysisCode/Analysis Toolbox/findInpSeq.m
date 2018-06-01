function [alphabet] = findInpSeq(asdf, delays, adj_mat, tol)

[~,N] = size(delays);
alphabet = cell(N,1);
parfor ii=1:N
    inpInds = find(adj_mat(:, ii));
    post_syn = asdf{ii};
    word_bag = cell(length(post_syn),1);
    if isempty(post_syn)
        continue;
    end
    for jj=1:length(inpInds)
        kk=1;
        inpI = inpInds(jj);
        dly = delays(inpI, ii) * asdf{end-1}(1);
        pre_syn = asdf{inpI} + dly;
        if isempty(pre_syn)
            continue;
        end
        for ll=1:length(post_syn)
            t_diff = (pre_syn(kk)-post_syn(ll));
            while kk<(length(pre_syn)-1) && t_diff <= 0
                if -t_diff < tol
                    sizeOrig = size(word_bag{ll},1);
                    word_bag{ll}(1+sizeOrig,1) = inpI;
                    word_bag{ll}(1+sizeOrig,2) =  ...
                        post_syn(ll) - (pre_syn(kk)-dly);
                    if word_bag{ll}(size(word_bag{ll},1),1) == 0
                        disp('huh');
                    end
                end
                kk=kk+1;
                t_diff = (pre_syn(kk)-post_syn(ll));
            end
        end
    end
    for jj=1:length(post_syn)
        if ~isempty(word_bag{jj})
%             if sum(word_bag{jj}(:,1)== 0)~= 0
%                 disp('huh');
%             end
            [~, I] = sort(word_bag{jj}(:,2), 'ascend');
            word_bag{jj} = word_bag{jj}(I,:);
        end
    end
    alphabet{ii} = word_bag;
    ii
end



end