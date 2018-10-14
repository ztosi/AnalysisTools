function [jsdsM, msesM] = distSimilarityTest(asdf)
isis = cellfun(@diff, asdf(1:end-2), 'UniformOutput', 0);
%isiMana = cellfun(@(x) diff(x(x>(max(x)/2))), asdf(1:end-2), 'UniformOutput', 0);
%randSelect = randperm(length(isiMana), 300);
%isiMana = isiMana(randSelect);
mxmx = max(cellfun(@max, isis))
mimi = min(cellfun(@min, isis))

disp(['Range: ' num2str(mxmx-mimi)]);


msesM = zeros(length(isis));
jsdsM = zeros(length(isis));


for ii=1:length(isis)
    parfor jj=(ii+1):length(isis)
        [jsdsM(ii,jj), msesM(ii,jj)] =  findMaxDistSimilarity2(isis{ii}, isis{jj}, ...
            'FullRange', [mimi mxmx], 'Transform', 'log10');
    end
    ii
end


msesM = msesM + msesM';
jsdsM = jsdsM + jsdsM';