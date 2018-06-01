IDX1 = kmeans(ydataNorm,7);
[idxSrt1, kMeansGrouping1 ]= sort(IDX1, 'descend');

figure; 
subplot(2,3,1);
scatter3(ydataNorm(:,1), ydataNorm(:,2), ydataNorm(:,3), 20, IDX1, 'filled');
title({'Normalized ssJSD', 'Color is  k-Means (k=7)  (Audio 2)'});
view([-101 25]);

subplot(2,3,2);
scatter3(ydataNorm(:,1), ydataNorm(:,2), ydataNorm(:,3), 20, modCol, 'filled');
title({'Normalized ssJSD', 'Color is Module (Audio 2)'});
view([-101 25]);

subplot(2,3,3);
s=scatter3(ydataNorm(:,1), ydataNorm(:,2), ydataNorm(:,3), 20, ei_c, 'filled');
%s.MarkerFaceAlpha=0.5;
title({'Normalized ssJSD', 'Color is E/I = (Audio 2)'});
view([-101 25]);


IDX2 = kmeans(ydataNoNorm,7);
[idxSrt2, kMeansGrouping2 ]= sort(IDX2, 'descend');

subplot(2,3,4);
scatter3(ydataNoNorm (:,1), ydataNoNorm(:,2), ydataNoNorm (:,3), 20, IDX2, 'filled');
title({'ssJSD', 'Color is  k-Means (k=7)  (Audio 2)'});
view([84 15]);

subplot(2,3,5);
scatter3(ydataNoNorm (:,1), ydataNoNorm (:,2), ydataNoNorm(:,3), 20, modCol, 'filled');
title({'ssJSD', 'Color is Module  (Audio 2)'});
view([84 15]);

subplot(2,3,6);
scatter3(ydataNoNorm (:,1), ydataNoNorm (:,2), ydataNoNorm (:,3), 20, ei_c, 'filled');
title({'ssJSD', 'Color is E/I  (Audio 2)'});
view([84 15]);


figure; 

subplot(211);
colormap(flipud(gray)); imagesc(rast_c(1000000:1002000,kMeansGrouping1)'); hold on;
for ii=1:6
bords = find(diff(idxSrt1));
plot([0 2000], [bords(ii) bords(ii)], 'Color', [0.45 0.6 1.0]);
end
title('Raster Seg. by Normalized ssJSD k-means (Audio 2)');

subplot(212);
colormap(flipud(gray)); imagesc(rast_c(1000000:1002000,kMeansGrouping2)'); hold on;
for ii=1:6
bords = find(diff(idxSrt2));
plot([0 2000], [bords(ii) bords(ii)], 'Color', [0.45 0.6 1.0]);
end
title('Raster Seg. by Regular ssJSD k-means  (Audio 2)');