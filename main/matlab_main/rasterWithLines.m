function rasterWithLines(rast, idx)


    [idxSrt, ord]  = sort(idx);
     
    chInds = find(diff(idxSrt)~=0);
    
   % figure; 
    hold on;
   % colormap(flipud(bone));
    imagesc(rast(:, ord)');
    
    [N,M] = size(rast);
    
    for ii=1:length(chInds)
       plot3([1 N], [chInds(ii) chInds(ii)], [5 5], 'LineWidth', 2, ...
           'Color', [.1, .7, .3]);  
    end

    xlim([1 N]);
    ylim([1 M]);
    
    view(2);
    hold off;

end